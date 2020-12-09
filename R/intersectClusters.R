#' Intersect pre-defined clusters
#'
#' Intersect pre-defined clusters from multiple modalities, 
#' pruning out combinations that are poorly separated based on the within-cluster sum of squares (WCSS).
#'
#' @param clusters A list of factors or vectors of the same length.
#' Each element corresponds to one modality and contains the cluster assignments for the same set of cells.
#' @param coords A list of matrices of length equal to \code{clusters}.
#' Each element should have number of rows equal to the number of cells (e.g., a matrix of PC coordinates);
#' we generally expect this to have been used to generate the corresponding entry of \code{clusters}.
#' @param scale Numeric scalar specifying the scaling factor to apply to the limit on the WCSS for each modality.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how parallelization should be performed.
#'
#' @details
#' We intersect clusters by only considering two cells to be in the same \dQuote{output} cluster if they are also clustered together in each modality.
#' In other words, all cells with a particular combination of identities in \code{clusters} are assigned to a separate output cluster.
#'
#' The simplest implementation of the above idea suffers from noise in the cluster definitions that introduces combinations with very few cells.
#' We eliminate these by greedily merging pairs of combinations, starting with the pairs that minimize the gain in the WCSS.
#' In this process, we only consider pairs of combinations that share at least cluster across all modalities (to avoid merges across unrelated clusters).
#'
#' A natural stopping point for this merging process is when the WCSS of the output clustering exceeds the WCSS of the original clustering for any modality.
#' This aims to preserve the original clustering in each modality by preventing overly aggressive merges that would greatly increase the WCSS,
#' while reducing the complexity of the output clustering by ensuring that the variance explained is comparable.
#'
#' Users can increase the aggressiveness of the merging procedure by increasing \code{scale}, e.g., to 1.05 or 1.
#' This will scale up the limit on the WCSS, allowing more merges to be performed before termination.
#'
#' @return An integer vector of length equal to the number of cells, containing the assignments to the output clusters.
#' @author Aaron Lun
#'
#' @examples
#' mat1 <- matrix(rnorm(10000), ncol=20)
#' chosen <- 1:250
#' mat1[chosen,1] <- mat1[chosen,1] + 10
#' clusters1 <- kmeans(mat1, 5)$cluster
#' table(clusters1, chosen=mat1[,1] > 5)
#'
#' # Pretending we have some other data for the same cells, e.g., ADT.
#' mat2 <- matrix(rnorm(10000), ncol=20)
#' chosen <- c(1:125, 251:375)
#' mat2[chosen,2] <- mat2[chosen,2] + 10
#' clusters2 <- kmeans(mat2, 5)$cluster
#' table(clusters2, mat2[,2] > 5)
#' 
#' # Intersecting the clusters:
#' clusters3 <- intersectClusters(list(clusters1, clusters2), list(mat1, mat2))
#' table(clusters3, mat1[,1] > 5)
#' table(clusters3, mat2[,2] > 5)
#' 
#' @export
#' @importFrom Matrix sparseMatrix
#' @importFrom S4Vectors DataFrame selfmatch
#' @importFrom DelayedArray getAutoBPPARAM setAutoBPPARAM
#' @importFrom scuttle .bpNotSharedOrUp
#' @importFrom BiocParallel bpstart bpstop SerialParam
intersectClusters <- function(clusters, coords, scale=1, BPPARAM=SerialParam()) {
    old <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(old), add=TRUE)
    if (!.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    if (is.null(names(clusters))) {
        names(clusters) <- sprintf("mode%i", seq_along(clusters))
    }
    names(coords) <- names(clusters)
    nmodes <- length(clusters)

    all.clusters <- do.call(DataFrame, c(clusters, list(check.names=FALSE)))
    combinations <- as.character(selfmatch(all.clusters))

    # Computing the WCSS for each modality (for the original clustering and for
    # each cluster combination).
    limits <- numeric(nmodes)
    all.combo.var <- vector("list", nmodes)
    for (m in seq_along(limits)) {
        limits[m] <- sum(.compute_wcss(clusters[[m]], coords[[m]])) * scale
        all.combo.var[[m]] <- .compute_wcss(combinations, coords[[m]])
    }

    # Computing statistics all possible pairwise merges.
    keep <- !duplicated(combinations)
    ucombinations <- combinations[keep]
    uclusters <- all.clusters[keep,,drop=FALSE]
    indices <- gains <- vector("list", nrow(uclusters))

    for (u in seq_len(nrow(uclusters))) {
        cur.combo <- ucombinations[u]
        cur.clusters <- uclusters[u,,drop=FALSE]

        everything <- integer(nrow(uclusters))
        for (m in seq_len(nmodes)) {
            everything <- everything + as.integer(uclusters[,m] == cur.clusters[,m])
        }

        allowed.friends <- ucombinations[everything > 0 & seq_along(everything) < u]
        out <- .compute_wcss_delta(cur.cluster=cur.combo, allowed=allowed.friends, 
            combinations=combinations, coords=coords, combo.var=all.combo.var)
        indices[[u]] <- out$indices
        gains[[u]] <- out$gains
    }

    # Iterating until any of the target variances are exceeded.
    indices <- do.call(rbind, indices)
    gains <- do.call(rbind, gains)
    running <- vapply(all.combo.var, sum, 0)
    summarized <- .summarize_wcss_gain(gains, running, limits)

    while (nrow(indices)) {
        # Picking the first merge with the smallest WCSS gain _and_ fits uder the limit.
        found <- FALSE
        for (chosen in order(summarized)) {
            running2  <- mapply("+", running, unlist(gains[chosen,], use.names=FALSE), USE.NAMES=FALSE)
            if (all(running2 < limits)) {
                running <- running2
                found <- TRUE
                break
            }
        }
        if (!found) {
            break
        }

        cur.index <- indices[chosen,]
        new.cluster <- paste0(cur.index[,1], ",", cur.index[,2])
        in.new.cluster <- combinations == cur.index[,1] | combinations == cur.index[,2]
        combinations[in.new.cluster] <- new.cluster

        # Identifying all friends for either of the two constituent clusters.
        self <- unlist(cur.index, use.names=FALSE)
        has.one <- indices[,1] %in% self | indices[,2] %in% self
        allowed.friends <- unlist(indices[has.one,], use.names=FALSE)
        allowed.friends <- setdiff(unique(allowed.friends), self)

        # Wiping out the entries for the merged clusters. 
        indices <- indices[!has.one,]
        gains <- gains[!has.one,]
        summarized <- summarized[!has.one]

        # Updating the WCSS entries.
        for (m in seq_len(nmodes)) {
            mode.combo.var <- all.combo.var[[m]] 
            other.clusters <- setdiff(names(mode.combo.var), self)
            mode.combo.var <- mode.combo.var[other.clusters]

            new.combo.var <- .compute_wcss_solo(coords[[m]][in.new.cluster,,drop=FALSE])
            names(new.combo.var) <- new.cluster
            all.combo.var[[m]] <- c(mode.combo.var, new.combo.var)
        }

        # Comparing the new cluster to all remaining clusters.
        out <- .compute_wcss_delta(cur.cluster=new.cluster, allowed=allowed.friends, 
            combinations=combinations, coords=coords, combo.var=all.combo.var)
        indices <- rbind(indices, out$indices)
        gains <- rbind(gains, out$gains)
        summarized <- c(summarized, .summarize_wcss_gain(out$gains, running, limits))
    }

    as.integer(factor(combinations))
}

#' @importFrom S4Vectors splitAsList
.compute_wcss <- function(mclusters, mcoords) {
    by.cluster <- splitAsList(seq_along(mclusters), mclusters)
    variances <- numeric(length(by.cluster))
    names(variances) <- names(by.cluster)

    for (bx in seq_along(by.cluster)) {
        current <- by.cluster[[bx]]
        tmp <- mcoords[current,,drop=FALSE]
        variances[bx] <- .compute_wcss_solo(tmp)
    }

    variances
}

#' @importFrom DelayedMatrixStats colVars
.compute_wcss_solo <- function(coords) {
    if (nrow(coords) >= 2) {
        sum(colVars(coords)) * (nrow(coords) - 1)
    } else {
        0
    }
}

#' @importFrom S4Vectors DataFrame 
.compute_wcss_delta <- function(cur.cluster, allowed, combinations, coords, combo.var) {
    gains <- vector("list", length(coords))
    names(gains) <- names(coords)

    for (m in seq_along(gains)) {
        mode.coords <- coords[[m]]
        mode.combo.var <- combo.var[[m]]
        gains[[m]] <- vapply(allowed, FUN=function(f) {
            to.merge <- c(cur.cluster, f)
            keep <- combinations %in% to.merge
            wcss <- .compute_wcss_solo(mode.coords[keep,,drop=FALSE])
            wcss - sum(mode.combo.var[to.merge])
        }, FUN.VALUE=0)
    }

    list(
        indices=DataFrame(first=rep(cur.cluster, length(allowed)), second=allowed),
        gains=do.call(DataFrame, gains) 
    )
}

.summarize_wcss_gain <- function(gains, running, limits) {
    # Normalizing each gain by the amount of remaining space in each mode,
    # and then taking the maximum normalized gain; by minimizing the maximum,
    # we aim to take the smallest step possible across all modes.
    remaining <- limits - running
    out <- mapply(`/`, as.list(gains), remaining, SIMPLIFY=FALSE)
    do.call(pmax, out)
}
