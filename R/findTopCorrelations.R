#' Find top correlations between genes
#'
#' For each gene, find the subset of other genes that have strongest positive/negative Spearman's rank correlations in a normalized expression matrix.
#'
#' @param x A numeric matrix containing a (possibly log-transformed) normalized expression matrix for genes (rows) and cells (columns).
#' Alternatively, a \linkS4class{SummarizedExperiment} containing such a matrix.
#' @param number Integer scalar specifying the number of top correlated genes to report for each gene in \code{x}.
#' @param y Another normalized expression matrix with the same number of cells as \code{x}, possibly with different features.
#' For the SummarizedExperiment method, this may also be a SummarizedExperiment object containing such a matrix.
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
#' Each DataFrame has up to \code{nrow(x) * number} rows, containing the top \code{number} correlated genes for each gene in \code{x}.
#' This contains the following fields:
#' \itemize{
#' \item \code{gene1}, the name (character) or row index (integer) of each gene in \code{x}.
#' \item \code{gene2}, the name (character) or row index (integer) of one of the top correlated genes to \code{gene1}.
#' This is another gene in \code{x} if \code{y=NULL}, otherwise it is a gene in \code{y}.
#' \item \code{rho}, the Spearman rank correlation for the current pair of \code{gene1} and \code{gene2}.
#' \item \code{p.value}, the approximate p-value associated with \code{rho} under the null hypothesis that the correlation is zero.
#' \item \code{FDR}, the adjusted p-value.
#' }
#' The rows are sorted by \code{gene1} and then \code{p.value}.
#' 
#' @details
#' In most cases, we only care about the top-correlated genes, allowing us to skip a lot of unnecessary computation.
#' This is achieved by transforming the problem of finding the largest Spearman correlation into a nearest-neighbor search in rank space.
#' For the sake of speed, we approximate the search by performing PCA to compress the rank values for all genes.
#'
#' We compute the p-value for each gene using the approximate method implemented in \code{\link{cor.test}}.
#' The FDR correction is performed by considering all possible pairs of genes, as these are implicitly tested in the neighbor search.
#' Note that this is somewhat conservative as it does not consider strong correlations outside the reported genes.
#'
#' @author Aaron Lun
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#'
#' df <- findTopCorrelations(sce, number=20) # top 20 correlated genes for each gene.
#' 
#' 
#' @name findTopCorrelations
NULL

######################
##### Internals ######
######################

#' @importFrom BiocSingular IrlbaParam
#' @importFrom BiocNeighbors findKNN queryKNN KmknnParam
#' @importFrom BiocParallel SerialParam bpstart bpstop
#' @importFrom scuttle .bpNotSharedOrUp
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

    if (is.null(y)) {
        .find_self_top_correlations(x, number=number, d=d, deferred=deferred, equiweight=equiweight,
            block=block, direction=direction, use.names=use.names,
            BSPARAM=BSPARAM, BNPARAM=BNPARAM, BPPARAM=BPPARAM)

    } else {
        .find_cross_top_correlations(x, y=y, number=number, d=d, deferred=deferred, equiweight=equiweight,
            block=block, direction=direction, use.names=use.names,
            BSPARAM=BSPARAM, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
    }
}

#' @importFrom Matrix t
#' @importFrom BiocNeighbors findKNN queryKNN buildIndex
#' @importFrom DelayedArray is_sparse
#' @importFrom beachmat rowBlockApply
#' @importFrom DelayedMatrixStats rowAnys
#' @importFrom S4Vectors List
.find_self_top_correlations <- function(x, direction, number, d, block, deferred, equiweight, use.names, BSPARAM, BNPARAM, BPPARAM) {
    if (direction=="positive") {
        combined <- t(x)
    } else {
        combined <- t(x)
        combined <- cbind(combined, -combined)
    }

    stash <- .create_blocked_rank_matrix(combined, deferred=deferred, block=block, d=d, 
        equiweight=equiweight, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
    rank.out <- stash$rank.out
    search.out <- stash$search.out
    nblocks <- stash$nblocks

    if (direction=="positive") {
        target <- search.out
        rank.x <- other.x <- rank.out
    } else {
        in.first <- seq_len(nrow(x))
        in.second <- nrow(x) + seq_len(nrow(x))
        target <- search.out[in.first,,drop=FALSE]
        other.x <- rank.out[in.second,,drop=FALSE]
        rank.x <- rank.out[in.first,,drop=FALSE]
    }
    precomputed <- buildIndex(target, BNPARAM=BNPARAM)

    output <- List()

    args <- list(
        ntests=nrow(x) * (nrow(x) - 1L), # don't divide by 2, as some pairs are tested twice if they show up reciprocally.
        names1=rownames(x),
        names2=rownames(x),
        use.names=use.names,
        nblocks=nblocks,
        equiweight=equiweight
    )

    if (direction %in% c("positive", "both")) {
        nn.out <- findKNN(BNINDEX=precomputed, k=number, BPPARAM=BPPARAM, get.distance=FALSE)
        rho <- rowBlockApply(rank.x, FUN=.compute_exact_neighbor_rho, other=rank.x, 
            nblocks=nblocks, indices=nn.out$index, BPPARAM=BPPARAM)
        output$positive <- do.call(.create_output_dataframe, c(list(nn.out$index, rho, positive=TRUE), args))
    }

    if (direction %in% c("negative", "both")) {
        # Plus 1 and then stripping out self, if it should come to that.
        nn.out <- queryKNN(query=search.out[in.second,,drop=FALSE],
            BNINDEX=precomputed, k=number+1L, get.distance=FALSE, BPPARAM=BPPARAM)

        discard <- nn.out$index==seq_len(nrow(nn.out$index))
        discard[!rowAnys(discard),ncol(discard)] <- TRUE # removing the last.
        stripped <- matrix(t(nn.out$index)[t(!discard)], ncol(nn.out$index) - 1L, nrow(nn.out$index)) # transposing so we refill by column-major.
        nn.out$index <- t(stripped)

        rho <- rowBlockApply(rank.x, FUN=.compute_exact_neighbor_rho, other=rank.x, 
            nblocks=nblocks, indices=nn.out$index, BPPARAM=BPPARAM)
        output$negative <- do.call(.create_output_dataframe, c(list(nn.out$index, rho, positive=FALSE), args))
    }

    output
}

#' @importFrom Matrix t
#' @importFrom BiocNeighbors queryKNN buildIndex
#' @importFrom DelayedArray is_sparse
#' @importFrom beachmat rowBlockApply
#' @importFrom S4Vectors List
.find_cross_top_correlations <- function(x, y, direction, number, d, block, deferred, equiweight, use.names, BSPARAM, BNPARAM, BPPARAM) {
    if (direction=="positive") {
        in.first <- seq_len(nrow(x))
        in.second <- nrow(x) + seq_len(nrow(y))
        combined <- cbind(t(x), t(y))
    } else {
        in.first <- seq_len(nrow(x))
        in.first.neg <- nrow(x) + in.first
        in.second <- 2L*nrow(x) + seq_len(nrow(y))

        # We make positive and negative versions of all of them for symmetry's sake.
        # Technically we only need a negative version of the smaller one... oh well.
        combined <- t(x)
        alt <- t(y)
        combined <- cbind(combined, -combined, alt, -alt)
    }

    stash <- .create_blocked_rank_matrix(combined, deferred=deferred, block=block, d=d, 
        equiweight=equiweight, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
    rank.out <- stash$rank.out
    search.out <- stash$search.out
    nblocks <- stash$nblocks

    precomputed <- buildIndex(search.out[in.second,,drop=FALSE], BNPARAM=BNPARAM)

    output <- List()

    args <- list(
        ntests=nrow(x) * nrow(y), 
        names1=rownames(x),
        names2=rownames(y),
        use.names=use.names,
        nblocks=nblocks,
        equiweight=equiweight
    )

    if (direction %in% c("positive", "both")) {
        nn.out <- queryKNN(query=search.out[in.first,,drop=FALSE], k=number, 
            get.distance=FALSE, BNINDEX=precomputed, BPPARAM=BPPARAM)
        rho <- rowBlockApply(rank.out[in.first,,drop=FALSE], FUN=.compute_exact_neighbor_rho, 
            nblocks=nblocks, other=rank.out[in.second,,drop=FALSE], indices=nn.out$index, BPPARAM=BPPARAM)
        output$positive <- do.call(.create_output_dataframe, c(list(nn.out$index, rho, positive=TRUE), args))
    }

    if (direction %in% c("negative", "both")) {
        nn.out <- queryKNN(query=search.out[in.first.neg,,drop=FALSE], k=number, 
            get.distance=FALSE, BNINDEX=precomputed, BPPARAM=BPPARAM)
        rho <- rowBlockApply(rank.out[in.first,,drop=FALSE], FUN=.compute_exact_neighbor_rho, 
            nblocks=nblocks, other=rank.out[in.second,,drop=FALSE], indices=nn.out$index, BPPARAM=BPPARAM)
        output$negative <- do.call(.create_output_dataframe, c(list(nn.out$index, rho, positive=FALSE), args))
    }

    output
}

#' @importFrom BiocGenerics cbind
#' @importFrom DelayedMatrixStats colVars
.create_blocked_rank_matrix <- function(x, deferred, block, ..., d, equiweight, BSPARAM, BPPARAM) {
    if (is.null(block)) {
        rank.out <- .create_rank_matrix(x, deferred=deferred, BPPARAM=BPPARAM)
        search.out <- .compress_rank_matrix(rank.out, d=d, deferred=deferred, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
        nblocks <- nrow(x)
    } else {
        # Remember, genes are columns coming into this function!
        # So we have to split by rows, as cells are in different batches.
        by.block <- split(seq_len(nrow(x)), block)
        rank.out <- search.out <- vector("list", length(by.block))

        for (i in seq_along(by.block)) {
            chosen <- by.block[[i]]
            rank.out[[i]] <- .create_rank_matrix(x[chosen,,drop=FALSE], deferred=deferred, BPPARAM=BPPARAM)
            search.out[[i]] <- .compress_rank_matrix(rank.out[[i]], d=d, deferred=deferred, BSPARAM=BSPARAM, BPPARAM=BPPARAM)

            # Equalizing the distance contribution based on the total variance.
            # This is not necessary for d=NA because the rank scaling does this for us,
            # but I do it again anyway because the PCA may change the relative contributions.
            if (equiweight) {
                total.var <- sum(colVars(search.out[[i]]))
                search.out[[i]] <- search.out[[i]] / sqrt(total.var)
            }
        }

        rank.out <- do.call(cbind, rank.out)
        nblocks <- lengths(by.block)
        search.out <- do.call(cbind, search.out)
    }

    list(rank.out=rank.out, search.out=search.out, nblocks=nblocks)
}

#' @importFrom Matrix colMeans t
#' @importFrom BiocSingular DeferredMatrix
#' @importFrom DelayedArray getAutoBPPARAM setAutoBPPARAM
#' @importFrom scran scaledColRanks
.create_rank_matrix <- function(x, deferred, ..., BPPARAM) {
    if (!deferred) {
        scaledColRanks(x, ..., transposed=TRUE)
    } else {
        old <- getAutoBPPARAM()
        setAutoBPPARAM(BPPARAM)
        on.exit(setAutoBPPARAM(old))

        y <- scaledColRanks(x, ..., transposed=FALSE, as.sparse=TRUE, BPPARAM=BPPARAM)
        y <- DeferredMatrix(y, center=colMeans(y))
        t(y)
    }
}

#' @importFrom BiocSingular runPCA
.compress_rank_matrix <- function(rank.x, d, deferred, BSPARAM, BPPARAM) {
    if (!is.na(d)) {
        BSPARAM@deferred <- deferred # TODO: add easier setters.
        runPCA(rank.x, rank=d, get.rotation=FALSE, BSPARAM=BSPARAM, BPPARAM=BPPARAM)$x
    } else {
        as.matrix(rank.x)
    }
}

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

#' @importFrom stats p.adjust
#' @importFrom metapod parallelStouffer
#' @importFrom S4Vectors DataFrame
.create_output_dataframe <- function(indices, rho, nblocks, equiweight, use.names, names1, names2, ntests, positive) {
    rho <- do.call(mapply, c(list(FUN=rbind, SIMPLIFY=FALSE), rho)) # still fragmented from the blockApply.

    if (length(nblocks)==1L) {
        mean.rho <- rho[[1]]
    } else if (equiweight) {
        mean.rho <- Reduce("+", rho)/length(rho)
    } else {
        mean.rho <- mapply("*", rho, nblocks, SIMPLIFY=FALSE)
        mean.rho <- Reduce("+", mean.rho)/sum(nblocks)
    }

    self <- rep(seq_len(nrow(indices)), ncol(indices))
    df <- DataFrame(gene1=self, gene2=as.vector(indices), rho=as.vector(mean.rho))

    # Assembling the DF.
    if (use.names) {
        if (!is.null(names1)) {
            df$gene1 <- names1[df$gene1]
        }
        if (!is.null(names2)) {
            df$gene2 <- names2[df$gene2]
        } 
    }

    # Mildly adapted from cor.test.
    p.values <- vector("list", length(rho))
    for (i in seq_along(rho)) {
        ncells <- nblocks[i]
        cur.rho <- as.vector(rho[[i]])

        q <- (ncells^3 - ncells) * (1 - cur.rho)/6
        den <- (ncells * (ncells^2 - 1)/6)
        r <- 1 - q/den
        tstat <- r/sqrt((1 - r^2)/(ncells - 2))

        p.values[[i]] <- pt(tstat, df = ncells - 2, lower.tail = !positive)
    }

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

    df$p.value <- p.value
    df$FDR <- p.adjust(df$p.value, method="BH", n = ntests)

    o <- order(df$gene1, df$p.value)
    df[o,,drop=FALSE]
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
