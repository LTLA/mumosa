#' Rescale matrices for different modes
#'
#' Rescale matrices for different data modalities so that their distances are more comparable, using the distances to neighbors to approximate noise.
#'
#' @param x A list of numeric matrices where each row is a cell and each column is some dimension/variable.
#' For gene expression data, this is usually the matrix of PC coordinates.
#' All matrices should have the same number of rows.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} containing relevant matrices in its assays.
#'
#' Alternatively, a \linkS4class{SingleCellExperiment} containing relevant matrices in its assays, \code{\link{reducedDims}} or \code{\link{altExps}}.
#' @param weights A numeric vector of length equal to \code{x} (if a list), specifying the weight of each mode.
#' Defaults to equal weights for all modes.
#' See details for how to interpret this argument when \code{x} is a SummarizedExperiment.
#' @param combine A logical scalar specifying whether the rescaled matrices should be combined into a single matrix. 
#' @inheritParams runMultiUMAP
#' @param assays A character or integer vector of assays to extract and transpose for use in the ANY method.
#' For the SingleCellExperiment, this argument can be missing, in which case no assays are used.
#' @param dimreds A character or integer vector of \code{\link{reducedDims}} to extract for use in the ANY method.
#' This argument can be missing, in which case no assays are used.
#' @param altexps A character or integer vector of \code{\link{altExps}} to extract and transpose for use in the ANY method.
#' This argument can be missing, in which case no alternative experiments are used.
#' @param ... For the generic, further arguments to pass to specific methods.
#'
#' For the SummarizedExperiment and SingleCellExperiment methods, further arguments to pass to the ANY method.
#' @param k An integer scalar specifying the number of neighbors to use for the distance calculation.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the algorithm to use for the nearest-neighbor search.
#' @param num.threads Integer scalar specifying the number of threads to use for the neighbor search.
#' @param BPPARAM Deprecated, use \code{num.threads} instead.
#' 
#' @return A numeric matrix with number of rows equal to the number of cells,
#' where the columns span all variables across all modes supplied in \code{x}.
#' Values are scaled so that each mode contributes the specified weight to downstream Euclidean distance calculations.
#'
#' If \code{combine=FALSE}, a list of rescaled matrices is returned instead.
#'
#' @details
#' When dealing with multi-modal data, we may wish to combine all modes into a single matrix for downstream processing.
#' However, a naive \code{cbind} does not account for the fact that different modes may very different scales and number of features.
#' A mode with a larger scale or more features may dominate steps such as clustering or dimensionality reduction.
#' This function attempts to rescale the contents for each matrix so that the modes are more comparable. 
#'
#' A naive approach to rescaling would be to just equalize the total variances across matrices.
#' This is not ideal as it fails to consider the differences in biological variation captured by each mode.
#' For example, if a biological phenomenon is only present in one mode, that matrix's total variance would naturally be higher.
#' Scaling all matrices to the same total variance would suppress genuine variation and inflate the relative contribution of noise.
#'
#' We instead use the distance to the \code{k}th nearest neighbor as an estimate of the per-mode \dQuote{noise}.
#' Modes with more features or higher technical noise will have larger distances, 
#' and downscaling each matrix by the median distance will correct for differences between modes.
#' At the same time, by only considering the nearest neighbors, 
#' we avoid capturing (and inadvertently eliminating) variance due to mode-specific population structure.
#'
#' The default approach is to weight each mode equally during the rescaling process,
#' i.e., the median distance to the \code{k}th nearest neighbor will be equal for all modes after rescaling.
#' However, we can also set \code{weights} to control the fold-differences in the median distances.
#' For example, a weight of 2 for one mode would mean that its median distance after rescaling is twice as large as that from a mode with a weight of 1. 
#' This may be useful for prioritizing modes that are more likely to be important.
#' 
#' The correspondence between non-\code{NULL} \code{weights} and the modes is slightly tricky whe \code{x} is not a list.
#' If \code{x} is a SummarizedExperiment, the modes are ordered as: all entries in \code{assays} in the specified order, then all entries in \code{extras}.
#' If \code{x} is a SingleCellExperiment, the modes are ordered as: all entries in \code{assays} in the specified order, 
#' then all entries in \code{dimreds}, then all entries in \code{altexps}, and finally all entries in \code{extras}.
#'
#' @examples
#' # Mocking up a gene expression + ADT dataset:
#' library(scater)
#' exprs_sce <- mockSCE()
#' exprs_sce <- logNormCounts(exprs_sce)
#' exprs_sce <- runPCA(exprs_sce)
#' 
#' adt_sce <- mockSCE(ngenes=20) 
#' adt_sce <- logNormCounts(adt_sce)
#' altExp(exprs_sce, "ADT") <- adt_sce
#'
#' combined <- rescaleByNeighbors(exprs_sce, dimreds="PCA", altexps="ADT")
#' dim(combined)
#' 
#' @author Aaron Lun
#' @name rescaleByNeighbors
NULL

#' @importFrom stats median
#' @importFrom BiocNeighbors findDistance
#' @importFrom BiocParallel SerialParam
#' @importFrom scuttle .bpNotSharedOrUp
#' @importFrom BiocParallel bpstart bpstop
#' @importFrom DelayedArray DelayedArray
.rescale_modal_matrices <- function(x, k=50, weights=NULL, combine=TRUE, num.threads=1, BNPARAM=NULL, BPPARAM=NULL) {
    if (length(x)==0) {
        stop("'x' must contain one or more matrices")
    }

    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    if (is.null(weights)) {
        weights <- 1
    }
    weights <- rep(weights, length.out=length(x))

    for (i in seq_along(x)) {
        nn <- findDistance(as.matrix(x[[i]]), k=k, BNPARAM=BNPARAM, BPPARAM=BPPARAM, num.threads=num.threads)
        nn <- median(nn)/weights[i]
        x[[i]] <- x[[i]]/max(nn, 1e-8)
    }

    if (combine) {
        # Work around Bioconductor/DelayedArray#100 for the time being.
        if (any(vapply(x, is, class="DelayedArray", TRUE))) {
            x <- lapply(x, DelayedArray)
        }
        do.call(cbind, x)
    } else {
        x
    }
}

#' @export
#' @rdname rescaleByNeighbors
setGeneric("rescaleByNeighbors", function(x, ...) standardGeneric("rescaleByNeighbors"))

#' @export
#' @rdname rescaleByNeighbors
setMethod("rescaleByNeighbors", "ANY", .rescale_modal_matrices)

#' @export
#' @rdname rescaleByNeighbors
setMethod("rescaleByNeighbors", "SummarizedExperiment", function(x, assays, extras=list(), ...) {
    targets <- .collate_se_matrices(x, assays)
    .rescale_modal_matrices(c(targets, extras), ...)
})

#' @export
#' @rdname rescaleByNeighbors
setMethod("rescaleByNeighbors", "SingleCellExperiment", 
    function(x, assays=NULL, dimreds=NULL, altexps=NULL, altexp.assay="logcounts", extras=list(), ...) 
{
    targets <- .collate_sce_matrices(x, dimreds=dimreds, altexps=altexps, altexp.assay=altexp.assay)
    callNextMethod(x, assays=assays, extras=c(targets, extras), ...)
})
