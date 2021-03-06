% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findTopCorrelations.R
\name{findTopCorrelations}
\alias{findTopCorrelations}
\alias{findTopCorrelations,ANY-method}
\alias{findTopCorrelations,SummarizedExperiment-method}
\title{Find top correlations between features}
\usage{
findTopCorrelations(x, number, ...)

\S4method{findTopCorrelations}{ANY}(
  x,
  number = 10,
  y = NULL,
  d = 50,
  direction = c("both", "positive", "negative"),
  subset.cols = NULL,
  block = NULL,
  equiweight = TRUE,
  use.names = TRUE,
  deferred = TRUE,
  BSPARAM = IrlbaParam(),
  BNPARAM = KmknnParam(),
  BPPARAM = SerialParam()
)

\S4method{findTopCorrelations}{SummarizedExperiment}(
  x,
  number,
  y = NULL,
  use.names = TRUE,
  ...,
  assay.type = "logcounts"
)
}
\arguments{
\item{x, y}{Normalized expression matrices containing features in the rows and cells in the columns.
Each matrix should have the same set of columns but a different set of features, usually corresponding to different modes for the same cells.

Alternatively, \linkS4class{SummarizedExperiment} objects containing such a matrix.

Finally, \code{y} may be \code{NULL}, in which correlations are computed between features in \code{x}.}

\item{number}{Integer scalar specifying the number of top correlated features to report for each feature in \code{x}.}

\item{...}{For the generic, further arguments to pass to specific methods.

For the SummarizedExperiment method, further arguments to pass to the ANY method.}

\item{d}{Integer scalar specifying the number of dimensions to use for the approximate search via PCA.
If \code{NA}, no approximation of the rank values is performed prior to the search.}

\item{direction}{String specifying the sign of the correlations to search for.}

\item{subset.cols}{Vector indicating the columns of \code{x} (and \code{y}) to retain for computing correlations.}

\item{block}{A vector or factor of length equal to the number of cells, specifying the block of origin for each cell.}

\item{equiweight}{Logical scalar indicating whether each block should be given equal weight, if \code{block} is specified.
If \code{FALSE}, each block is weighted by the number of cells.}

\item{use.names}{Logical scalar specifying whether row names of \code{x} and/or \code{y} should be reported in the output, if available.

For the SummarizedExperiment method, this may also be a string specifying the \code{\link{rowData}} column containing the names to use;
or a character vector of length 2, where the first and second entries specify the \code{\link{rowData}} columns containing the names in \code{x} and \code{y} respectively.
If either entry is \code{NA}, the existing row names for the corresponding object are used.
Note that this only has an effect on \code{y} if it is a SummarizedExperiment.}

\item{deferred}{Logical scalar indicating whether a fast deferred calculation should be used for the rank-based PCA.}

\item{BSPARAM}{A \linkS4class{BiocSingularParam} object specifying the algorithm to use for the PCA.}

\item{BNPARAM}{A \linkS4class{BiocNeighborParam} object specifying the algorithm to use for the neighbor search.}

\item{BPPARAM}{A \linkS4class{BiocParallelParam} object specifying the parallelization scheme to use.}

\item{assay.type}{String or integer scalar specifying the assay containing the matrix of interest in \code{x} (and \code{y}, if a SummarizedExperiment).}
}
\value{
A \linkS4class{List} containing one or two \linkS4class{DataFrame}s for results in each direction.
These are named \code{"positive"} and \code{"negative"}, and are generated according to \code{direction};
if \code{direction="both"}, both DataFrames will be present.

Each DataFrame has up to \code{nrow(x) * number} rows, containing the top \code{number} correlated features for each feature in \code{x}.
This contains the following fields:
\itemize{
\item \code{feature1}, the name (character) or row index (integer) of each feature in \code{x}.
Not all features may be reported here, see Details.
\item \code{feature2}, the name (character) or row index (integer) of one of the top correlated features to \code{feature1}.
This is another feature in \code{x} if \code{y=NULL}, otherwise it is a feature in \code{y}.
\item \code{rho}, the Spearman rank correlation for the current pair of \code{feature1} and \code{feature2}.
\item \code{p.value}, the approximate p-value associated with \code{rho} under the null hypothesis that the correlation is zero.
\item \code{FDR}, the adjusted p-value.
}
The rows are sorted by \code{feature1} and then \code{p.value}.
}
\description{
For each feature, find the subset of other features in the same or another modality that have strongest positive/negative Spearman's rank correlations in a pair of normalized expression matrices.
}
\details{
In most cases, we only care about the top-correlated features, allowing us to skip a lot of unnecessary computation.
This is achieved by transforming the problem of finding the largest Spearman correlation into a nearest-neighbor search in rank space.
For the sake of speed, we approximate the search by performing PCA to compress the rank values for all features.

For each direction, we compute the one-sided p-value for each feature using the approximate method implemented in \code{\link{cor.test}}.
The FDR correction is performed by considering all possible pairs of features, as these are implicitly tested in the neighbor search.
Note that this is somewhat conservative as it does not consider strong correlations outside the reported features.

If \code{block} is specified, correlations are computed separately for each block of cells.
For each feature pair, the reported \code{rho} is set to the average of the correlations across all blocks.
Similarly, the p-value corresponding to each correlation is computed separately for each block and then combined across blocks with Stouffer's method.
If \code{equiweight=FALSE}, the average correlation and each per-block p-value is weighted by the number of cells.

We only consider pairs of features that have computable correlations in at least one block.
Blocks are ignored if one or the other feature has tied values (typically zeros) for all cells in that block.
This means that a feature may not have any entries in \code{feature1} if it forms no valid pairs, e.g., because it is not expressed.
Similarly, the total number of rows may be less than the maximum if insufficient valid pairs are available.
}
\examples{
library(scuttle)
sce1 <- mockSCE()
sce1 <- logNormCounts(sce1)

sce2 <- mockSCE(ngenes=20) # pretend this is CITE-seq data, or something.
sce2 <- logNormCounts(sce2)

# Top 20 correlated features in 'sce2' for each feature in 'sce1':
df <- findTopCorrelations(sce1, sce2, number=20) 
df

}
\seealso{
\code{\link{computeCorrelations}}, to compute correlations for all pairs of features.
}
\author{
Aaron Lun
}
