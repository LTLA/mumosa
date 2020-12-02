#' Find top correlations between features
#'
#' For each feature, find the subset of other features in the same or another modality that have strongest positive/negative Spearman's rank correlations in a pair of normalized expression matrices.
#'
#' @param x,y Normalized expression matrices containing features in the rows and cells in the columns.
#' Each matrix should have the same set of columns but a different set of features, usually corresponding to different modes for the same cells.
#' 
#' Alternatively, \linkS4class{SummarizedExperiment} objects containing such a matrix.
#'
#' \code{y} may be missing, in which correlations are computed between features in \code{x}.
#' @param number Integer scalar specifying the number of top correlated features to report for each feature in \code{x}.
#' @param block A vector or factor of length equal to the number of cells, specifying the block of origin for each cell.
#' @param d Integer scalar specifying the number of dimensions to use for the approximate search via PCA.
#' If \code{NA}, no approximation of the rank values is performed prior to the search.
#' @param use.names Logical scalar specifying whether row names of \code{x} and/or \code{y} should be reported in the output, if available.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying the algorithm to use for the PCA.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the algorithm to use for the neighbor search.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying the parallelization scheme to use.
#' @param ... For the generic, further arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#' @param assay.type String or integer scalar specifying the assay containing the matrix of interest in \code{x} (and \code{y}, if a SummarizedExperiment).
#' @param direction String specifying the sign of the correlations to search for.
#'
#' @return A \linkS4class{List} containing one or two \linkS4class{DataFrame}s for results in each direction.
#' These are named \code{"positive"} and \code{"negative"}, and are generated according to \code{direction};
#' if \code{direction="both"}, both DataFrames will be present.
#'
#' Each DataFrame has up to \code{nrow(x) * number} rows, containing the top \code{number} correlated features for each feature in \code{x}.
#' This contains the following fields:
#' \itemize{
#' \item \code{feature1}, the name (character) or row index (integer) of each feature in \code{x}.
#' \item \code{feature2}, the name (character) or row index (integer) of one of the top correlated features to \code{feature1}.
#' This is another feature in \code{x} if \code{y=NULL}, otherwise it is a feature in \code{y}.
#' \item \code{rho}, the Spearman rank correlation for the current pair of \code{feature1} and \code{feature2}.
#' \item \code{p.value}, the approximate p-value associated with \code{rho} under the null hypothesis that the correlation is zero.
#' \item \code{FDR}, the adjusted p-value.
#' }
#' The rows are sorted by \code{feature1} and then \code{p.value}.
#' 
#' @details
#' In most cases, we only care about the top-correlated features, allowing us to skip a lot of unnecessary computation.
#' This is achieved by transforming the problem of finding the largest Spearman correlation into a nearest-neighbor search in rank space.
#' For the sake of speed, we approximate the search by performing PCA to compress the rank values for all features.
#'
#' For each direction, we compute the one-sided p-value for each feature using the approximate method implemented in \code{\link{cor.test}}.
#' The FDR correction is performed by considering all possible pairs of features, as these are implicitly tested in the neighbor search.
#' Note that this is somewhat conservative as it does not consider strong correlations outside the reported features.
#'
#' If \code{block} is specified, correlations are computed separately for each block of cells.
#' For each feature pair, the reported \code{rho} is set to the average of the correlations across all blocks.
#' Similarly, the p-value corresponding to each correlation is computed separately for each block and then combined across blocks with Stouffer's method.
#' If \code{equiweight=FALSE}, the average correlation and each per-block p-value is weighted by the number of cells.
#'
#' @author Aaron Lun
#'
#' @examples
#' library(scuttle)
#' sce1 <- mockSCE()
#' sce1 <- logNormCounts(sce1)
#'
#' sce2 <- mockSCE(ngenes=20) # pretend this is CITE-seq data, or something.
#' sce2 <- logNormCounts(sce2)
#'
#' # Top 20 correlated features in 'sce2' for each feature in 'sce1':
#' df <- findTopCorrelations(sce1, sce2, number=20) 
#' df
#' 
#' @seealso
#' \code{\link{computeCorrelations}}, to compute correlations for all pairs of features.
#' 
#' @name findTopCorrelations
NULL

#' @importFrom S4Vectors DataFrame splitAsList
#' @importFrom IRanges heads
#' @importFrom stats p.adjust
#' @importFrom BiocSingular IrlbaParam
#' @importFrom BiocNeighbors findKNN queryKNN KmknnParam buildIndex
#' @importFrom BiocParallel SerialParam bpstart bpstop
#' @importFrom scuttle .bpNotSharedOrUp
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats rowAnys
.find_top_correlations <- function(x, number=10, y=NULL, d=50, 
    direction=c("both", "positive", "negative"), 
    block=NULL, equiweight=TRUE, use.names=TRUE, deferred=TRUE, 
    BSPARAM=IrlbaParam(), BNPARAM=KmknnParam(), BPPARAM=SerialParam()) 
{
    direction <- match.arg(direction)
    if (!.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    self.search <- is.null(y)

    # Wrapping in DA's to avoid actually subsetting the matrices unnecessarily.
    x <- DelayedArray(x)
    if (self.search) {
        y <- x
    } else {
        y <- DelayedArray(y)
    }

    xcat <- .find_block_availability(x, block=block, BPPARAM=BPPARAM)
    if (self.search) {
        ycat <- xcat
    } else {
        ycat <- .find_block_availability(y, block=block, BPPARAM=BPPARAM)
    }

    total.ntests <- 0L
    output <- List()
    empty <- DataFrame(feature1=integer(0), feature2=integer(0), rho=numeric(0), p.value=numeric(0))
    if (direction %in% c("positive", "both")) {
        output$positive <- list(empty)
    }
    if (direction %in% c("negative", "both")) {
        output$negative <- list(empty)
    }
    ocounter <- 1L

    for (x.chosen in xcat$by) {
        x.valid <- unlist(xcat$available[x.chosen[1],,drop=FALSE])
        for (y.chosen in ycat$by) {
            y.valid <- unlist(ycat$available[y.chosen[1],,drop=FALSE])
            both.valid <- x.valid & y.valid
            if (!any(both.valid)) {
                next
            }

            same <- self.search && identical(x.chosen, y.chosen)
            if (same) {
                ntests <- length(x.chosen) * (length(x.chosen) - 1L) # don't divide by 2, as the same pair may be reported twice.
            } else {
                ntests <- length(x.chosen) * length(y.chosen)
            }
            total.ntests <- total.ntests + ntests

            block.info <- xcat$block[both.valid]
            cells <- unlist(block.info)
            x0 <- x[x.chosen,cells,drop=FALSE]
            if (!same) {
                y0 <- y[y.chosen,cells,drop=FALSE]
            }

            if (direction=="positive") {
                in.first <- seq_len(nrow(x0))
                if (same) {
                    in.second <- in.first
                    combined <- t(x0)
                } else {
                    in.second <- nrow(x0) + seq_len(nrow(y0))
                    combined <- cbind(t(x0), t(y0))
                }
            } else {
                in.first <- seq_len(nrow(x0))
                in.first.neg <- nrow(x0) + in.first
                combined <- t(x0)

                if (same) {
                    in.second <- in.first
                    combined <- cbind(combined, -combined)
                } else {
                    in.second <- 2L*nrow(x0) + seq_len(nrow(y0))

                    # We make positive and negative versions of all of them for symmetry's sake.
                    # Technically we only need a negative version of the smaller one... oh well.
                    alt <- t(y0)
                    combined <- cbind(combined, -combined, alt, -alt)
                }
            }

            nblocks <- lengths(block.info)
            stash <- .create_blocked_rank_matrix(combined, deferred=deferred, nblocks=nblocks,
                d=d, equiweight=equiweight, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
            rank.out <- stash$rank.out
            search.out <- stash$search.out

            precomputed <- buildIndex(search.out[in.second,,drop=FALSE], BNPARAM=BNPARAM)

            if (!is.null(output$positive)) {
                if (same) {
                    nn.out <- findKNN(BNINDEX=precomputed, k=number, BPPARAM=BPPARAM, get.distance=FALSE)
                } else {
                    nn.out <- queryKNN(query=search.out[in.first,,drop=FALSE], k=number, 
                        get.distance=FALSE, BNINDEX=precomputed, BPPARAM=BPPARAM)
                }

                rho <- rowBlockApply(rank.out[in.first,,drop=FALSE], FUN=.compute_exact_neighbor_rho, 
                    nblocks=nblocks, other=rank.out[in.second,,drop=FALSE], indices=nn.out$index, BPPARAM=BPPARAM)
                output$positive[[ocounter]] <- .create_output_dataframe(nn.out$index, rho, positive=TRUE, 
                    nblocks=nblocks, equiweight=equiweight, chosen1=x.chosen, chosen2=y.chosen)
            }

            if (!is.null(output$negative)) {
                nn.out <- queryKNN(query=search.out[in.first.neg,,drop=FALSE], k=number + as.integer(same), 
                    get.distance=FALSE, BNINDEX=precomputed, BPPARAM=BPPARAM)

                # Searching for number + 1 and then stripping out any self-matches.
                if (same) {
                    discard <- nn.out$index==seq_len(nrow(nn.out$index))
                    discard[!rowAnys(discard),ncol(discard)] <- TRUE # removing the last.
                    stripped <- matrix(t(nn.out$index)[t(!discard)], ncol(nn.out$index) - 1L, nrow(nn.out$index)) # transposing so we refill by column-major.
                    nn.out$index <- t(stripped)
                }

                rho <- rowBlockApply(rank.out[in.first,,drop=FALSE], FUN=.compute_exact_neighbor_rho, 
                    nblocks=nblocks, other=rank.out[in.second,,drop=FALSE], indices=nn.out$index, BPPARAM=BPPARAM)
                output$negative[[ocounter]] <- .create_output_dataframe(nn.out$index, rho, positive=FALSE, 
                    nblocks=nblocks, equiweight=equiweight, chosen1=x.chosen, chosen2=y.chosen)
            }

            ocounter <- ocounter + 1L
        }
    }

    for (dir in names(output)) {
        combined <- do.call(rbind, output[[dir]])
        combined <- combined[order(combined$feature1, combined$p.value),,drop=FALSE]
        fragged <- splitAsList(combined, combined$feature1)
        fragged <- heads(fragged, number)
        final <- unlist(fragged, use.names=FALSE)
        final <- .fill_names(final, use.names, rownames(x), rownames(y))
        final$FDR <- p.adjust(final$p.value, method="BH", n=total.ntests)
        output[[dir]] <- final
    }

    output
}

######################
##### Internals ######
######################

#' @importFrom S4Vectors DataFrame splitAsList
#' @importFrom BiocParallel SerialParam
#' @importFrom DelayedArray DelayedArray setAutoBPPARAM getAutoBPPARAM
#' @importFrom DelayedMatrixStats rowVars
.find_block_availability <- function(x, block, tol=1e-8, BPPARAM=SerialParam()) {
    old <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(old))

    x <- DelayedArray(x)
    if (is.null(block)) {
        output <- list(all=rowVars(x) > tol) 
        by.block <- list(seq_len(ncol(x)))
    } else {
        by.block <- splitAsList(seq_len(ncol(x)), block)
        output <- vector("list", length(by.block))
        names(output) <- names(by.block)
        for (i in seq_along(by.block)) {
            output[[i]] <- rowVars(x, cols=by.block[[i]]) > tol
        }
    }

    availabilities <- do.call(DataFrame, c(output, list(check.names=FALSE)))
    ids <- selfmatch(availabilities)
    by.availability <- split(seq_along(ids), ids)
    list(by=by.availability, block=by.block, available=availabilities)
}

#' @importFrom BiocGenerics cbind
#' @importFrom DelayedMatrixStats colVars
#' @importFrom Matrix colMeans t
#' @importFrom BiocSingular DeferredMatrix
#' @importFrom DelayedArray getAutoBPPARAM setAutoBPPARAM
#' @importFrom scran scaledColRanks
#' @importFrom BiocSingular runPCA
.create_blocked_rank_matrix <- function(x, deferred, nblocks, ..., d, equiweight, BSPARAM, BPPARAM) {
    old <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(old))

    # Remember, features are columns coming into this function!
    # So we have to split by rows, as cells are in different batches.
    rank.out <- search.out <- vector("list", length(nblocks))
    last <- 0L

    for (i in seq_along(nblocks)) {
        chosen <- last + seq_len(nblocks[i])
        x0 <- x[chosen,,drop=FALSE]

        if (!deferred) {
            rank.out[[i]] <- scaledColRanks(x0, BPPARAM=BPPARAM, transposed=TRUE)
        } else {
            y <- scaledColRanks(x0, transposed=FALSE, as.sparse=TRUE, BPPARAM=BPPARAM)
            y <- DeferredMatrix(y, center=colMeans(y))
            rank.out[[i]] <- t(y)
        }

        if (!is.na(d)) {
            BSPARAM@deferred <- deferred # TODO: add easier setters.
            search.out[[i]] <- runPCA(rank.out[[i]], rank=d, get.rotation=FALSE, BSPARAM=BSPARAM, BPPARAM=BPPARAM)$x
        } else {
            search.out[[i]] <- as.matrix(rank.out[[i]])
        }

        # Equalizing the distance contribution based on the total variance.
        # This is not necessary for d=NA because the rank scaling does this for us,
        # but the PCA may change the relative contributions, so we have to do something.
        if (!is.na(d) && equiweight) {
            total.var <- sum(colVars(search.out[[i]]))
            search.out[[i]] <- search.out[[i]] / sqrt(total.var)
        }

        last <- last + nblocks[i]
    }

    rank.out <- do.call(cbind, rank.out)
    search.out <- do.call(cbind, search.out)
    list(rank.out=rank.out, search.out=search.out)
}

#' @importFrom Matrix rowSums
#' @importFrom DelayedArray currentViewport makeNindexFromArrayViewport
.compute_exact_neighbor_rho <- function(block, other, indices, nblocks) {
    vp <- currentViewport()
    subset <- makeNindexFromArrayViewport(vp, expand.RangeNSBS=TRUE)
    if (!is.null(subset[[1]])) {
        indices <- indices[subset[[1]],,drop=FALSE]
    }

    output <- matrix(0, nrow(block), ncol(indices))
    output <- rep(list(output), length(nblocks))

    for (j in seq_len(ncol(indices))) {
        current <- block - as.matrix(other[indices[,j],,drop=FALSE])
        if (length(nblocks)==1L) {
            output[[1]][,j] <- 1 - 2*rowSums(current^2)
        } else {
            last <- 0L
            for (n in seq_along(nblocks)) {
                chosen <- last + seq_len(nblocks[n])
                output[[n]][,j] <- 1 - 2*rowSums(current[,chosen,drop=FALSE]^2)
                last <- last + nblocks[n]
            }
        }
    }

    output
}

#' @importFrom S4Vectors DataFrame
#' @importFrom metapod parallelStouffer
.create_output_dataframe <- function(indices, rho, nblocks, equiweight, positive, chosen1, chosen2) {
    rho <- do.call(mapply, c(list(FUN=rbind, SIMPLIFY=FALSE), rho)) # still fragmented from the blockApply.
    mean.rho <- .compute_mean_rho(rho, nblocks, equiweight)

    self <- rep(chosen1, ncol(indices))
    df <- DataFrame(feature1=self, feature2=chosen2[as.vector(indices)], rho=as.vector(mean.rho))

    p.values <- mapply(FUN=.compute_cor_p, rho=rho, ncells=nblocks, 
        MoreArgs=list(positive=positive), SIMPLIFY=FALSE)

    if (length(nblocks)==1L) {
        p.value <- p.values[[1]]
    } else {
        # Combining all the one-sided p-values together.
        if (equiweight) {
            weights <- NULL
        } else {
            weights <- nblocks
        }
        p.value <- parallelStouffer(p.values, weights=weights)$p.value
    }

    df$p.value <- as.vector(p.value)
    df
}

#######################
##### S4 methods ######
#######################

#' @export
#' @rdname findTopCorrelations
setGeneric("findTopCorrelations", function(x, number, ...) standardGeneric("findTopCorrelations"))

#' @export
#' @rdname findTopCorrelations
setMethod("findTopCorrelations", "ANY", .find_top_correlations)

#' @export
#' @rdname findTopCorrelations
#' @importFrom SummarizedExperiment assay
setMethod("findTopCorrelations", "SummarizedExperiment", function(x, number, y=NULL, ..., assay.type="logcounts") {
    if (is(y, "SummarizedExperiment")) {
        y <- assay(y, assay.type)
    }
    x <- assay(x, assay.type)
    .find_top_correlations(x, number=number, y=y, ...)
})
