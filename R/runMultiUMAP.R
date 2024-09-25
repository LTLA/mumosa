#' Multi-modal UMAP
#'
#' Perform UMAP with multiple input matrices by intersecting their simplicial sets.
#' Typically used to combine results from multiple data modalities into a single embedding.
#'
#' @param x For \code{calculateMultiUMAP}, a list of numeric matrices where each row is a cell and each column is some dimension/variable.
#' For gene expression data, this is usually the matrix of PC coordinates.
#' All matrices should have the same number of rows.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} containing relevant matrices in its assays.
#'
#' Alternatively, a \linkS4class{SingleCellExperiment} containing relevant matrices in its assays, \code{\link{reducedDims}} or \code{\link{altExps}}.
#' This is also the only permissible argument for \code{runMultiUMAP}.
#' @param assays A character or integer vector of assays to extract and transpose for use in the UMAP.
#' For the SingleCellExperiment, this argument can be missing, in which case no assays are used.
#' @param dimreds A character or integer vector of \code{\link{reducedDims}} to extract for use in the UMAP.
#' This argument can be missing, in which case no assays are used.
#' @param altexps A character or integer vector of \code{\link{altExps}} to extract and transpose for use in the UMAP.
#' This argument can be missing, in which case no alternative experiments are used.
#' @param altexp.assay A character or integer vector specifying the assay to extract from alternative experiments, when \code{altexp} is specified.
#' This is recycled to the same length as \code{altexp}.
#' @param extras A list of further matrices of similar structure to those matrices in a list-like \code{x}.
#' @param ... For the generic, further arguments to pass to specific methods.
#'
#' For the ANY method, further arguments to pass to \code{\link[uwot]{umap}}.
#'
#' For the SummarizedExperiment and SingleCellExperiment methods, and for \code{runMultiUMAP}, further arguments to pass to the ANY method.
#' @param metric Character vector specifying the type of distance to use for each matrix in \code{x}.
#' This is recycled to the same number of matrices supplied in \code{x}.
#' @param name String specifying the name of the \code{\link{reducedDims}} in which to store the UMAP.
#'
#' @return 
#' For \code{calculateMultiUMAP}, a numeric matrix containing the low-dimensional UMAP embedding.
#'
#' For \code{runMultiUMAP}, \code{x} is returned with a \code{MultiUMAP} field in its \code{\link{reducedDims}}.
#'
#' @details
#' These functions serve as convenience wrappers around \code{\link[uwot]{umap}} for multi-modal analysis.
#' The idea is that each input matrix in \code{x} corresponds to data for a different mode.
#' A typical example would consist of the PC coordinates generated from gene expression counts,
#' plus the log-abundance matrix for ADT counts from CITE-seq experiments;
#' one might also include matrices of transformed intensities from indexed FACS, to name some more possibilities.
#'
#' Roughly speaking, the idea is to identify nearest neighbors \emph{within} each mode to construct the simplicial sets.
#' Integration of multiple modes is performed by intersecting the sets to obtain a single graph, which is used in the rest of the UMAP algorithm.
#' By performing an intersection, we focus on relationships between cells that are consistently neighboring across all the modes,
#' thus providing greater resolution of differences at any mode.
#' The neighbor search within each mode also avoids difficulties with quantitative comparisons of distances between modes.
#' 
#' The most obvious use of this function is to generate a low-dimensional embedding for visualization.
#' However, users can also set \code{n_components} to a higher value (e.g., 10-20) to retain more information for downstream steps like clustering.
#' Do, however, remember to set the seed appropriately.
#'
#' By default, all modes use the distance metric of \code{metric} to construct the simplicial sets \emph{within} each mode.
#' However, it is possible to vary this by supplying a vector of metrics, e.g., \code{"euclidean"} for the first matrix, \code{"manhattan"} for the second.
#' For the SingleCellExperiment method, matrices are extracted in the order of assays, reduced dimensions and alternative experiments,
#' so any variation in \code{metrics} is also assumed to follow this order.
#'
#' @seealso
#' \code{\link[scater]{runUMAP}}, for the more straightforward application of UMAP.
#'
#' @author Aaron Lun
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
#' # Running a multimodal analysis using PCs for expression
#' # and log-counts for the ADTs. Annoyingly, have to prefix
#' # this for the time being to distinguish from the scater generic.
#' exprs_sce <- mumosa::runMultiUMAP(exprs_sce, dimreds="PCA", altexps="ADT")
#' plotReducedDim(exprs_sce, "MultiUMAP")
#' 
#' @name runMultiUMAP
NULL

#' @export 
#' @rdname runMultiUMAP
setGeneric("calculateMultiUMAP", function(x, ...) standardGeneric("calculateMultiUMAP"))

#' @export 
#' @rdname runMultiUMAP
#' @importFrom utils head
#' @importFrom uwot umap
#' @importFrom DelayedArray DelayedArray
setMethod("calculateMultiUMAP", "ANY", function(x, ..., metric="euclidean") {
    if (length(x)==0) {
        stop("'x' must contain one or more matrices")
    }
    mult.metrics <- .compute_multi_modal_metrics(x, metric=metric)

    # Work around Bioconductor/DelayedArray#100 for the time being.
    if (any(vapply(x, is, class="DelayedArray", TRUE))) {
        x <- lapply(x, DelayedArray)
    }

    combined <- as.matrix(do.call(cbind, x))
    umap(combined, metric=mult.metrics, ...)
})

.compute_multi_modal_metrics <- function(x, metric="euclidean") {
    ncols <- vapply(x, ncol, 0L)
    mult.metrics <- lapply(ncols, seq_len)
    mult.metrics <- mapply("+", mult.metrics, cumsum(c(0L, head(ncols, -1L))), SIMPLIFY=FALSE)
    names(mult.metrics) <- rep(metric, length.out=length(mult.metrics))
    mult.metrics
}

#' @export
#' @rdname runMultiUMAP
setMethod("calculateMultiUMAP", "SummarizedExperiment", function(x, assays, extras=list(), ...) {
    targets <- .collate_se_matrices(x, assays)
    callNextMethod(c(targets, extras), ...)
}) 

#' @export
#' @rdname runMultiUMAP
setMethod("calculateMultiUMAP", "SingleCellExperiment", 
    function(x, assays=NULL, dimreds=NULL, altexps=NULL, altexp.assay="logcounts", extras=list(), ...) 
{
    targets <- .collate_sce_matrices(x, dimreds=dimreds, altexps=altexps, altexp.assay=altexp.assay)
    callNextMethod(x, assays=assays, extras=c(targets, extras), ...)
}) 

#' @export
#' @rdname runMultiUMAP
#' @importFrom SingleCellExperiment reducedDim<-
runMultiUMAP <- function(x, ..., name="MultiUMAP") {
    reducedDim(x, name) <- calculateMultiUMAP(x, ...)
    x
}
