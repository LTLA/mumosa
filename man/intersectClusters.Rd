% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intersectClusters.R
\name{intersectClusters}
\alias{intersectClusters}
\title{Intersect pre-defined clusters}
\usage{
intersectClusters(clusters, coords, scale = 1, BPPARAM = SerialParam())
}
\arguments{
\item{clusters}{A list of factors or vectors of the same length.
Each element corresponds to one modality and contains the cluster assignments for the same set of cells.}

\item{coords}{A list of matrices of length equal to \code{clusters}.
Each element should have number of rows equal to the number of cells (e.g., a matrix of PC coordinates);
we generally expect this to have been used to generate the corresponding entry of \code{clusters}.}

\item{scale}{Numeric scalar specifying the scaling factor to apply to the limit on the WCSS for each modality.}

\item{BPPARAM}{A \linkS4class{BiocParallelParam} object specifying how parallelization should be performed.}
}
\value{
An integer vector of length equal to the number of cells, containing the assignments to the output clusters.
}
\description{
Intersect pre-defined clusters from multiple modalities, 
pruning out combinations that are poorly separated based on the within-cluster sum of squares (WCSS).
}
\details{
We intersect clusters by only considering two cells to be in the same \dQuote{output} cluster if they are also clustered together in each modality.
In other words, all cells with a particular combination of identities in \code{clusters} are assigned to a separate output cluster.

The simplest implementation of the above idea suffers from noise in the cluster definitions that introduces combinations with very few cells.
We eliminate these by greedily merging pairs of combinations, starting with the pairs that minimize the gain in the WCSS.
In this process, we only consider pairs of combinations that share at least cluster across all modalities (to avoid merges across unrelated clusters).

A natural stopping point for this merging process is when the WCSS of the output clustering exceeds the WCSS of the original clustering for any modality.
This aims to preserve the original clustering in each modality by preventing overly aggressive merges that would greatly increase the WCSS,
while reducing the complexity of the output clustering by ensuring that the variance explained is comparable.

Users can increase the aggressiveness of the merging procedure by increasing \code{scale}, e.g., to 1.05 or 1.
This will scale up the limit on the WCSS, allowing more merges to be performed before termination.
}
\examples{
mat1 <- matrix(rnorm(10000), ncol=20)
chosen <- 1:250
mat1[chosen,1] <- mat1[chosen,1] + 10
clusters1 <- kmeans(mat1, 5)$cluster
table(clusters1, chosen=mat1[,1] > 5)

# Pretending we have some other data for the same cells, e.g., ADT.
mat2 <- matrix(rnorm(10000), ncol=20)
chosen <- c(1:125, 251:375)
mat2[chosen,2] <- mat2[chosen,2] + 10
clusters2 <- kmeans(mat2, 5)$cluster
table(clusters2, mat2[,2] > 5)

# Intersecting the clusters:
clusters3 <- intersectClusters(list(clusters1, clusters2), list(mat1, mat2))
table(clusters3, mat1[,1] > 5)
table(clusters3, mat2[,2] > 5)

}
\author{
Aaron Lun
}
