# This script tests the rescaleByNeighbors functionality.
# library(testthat); library(mumosa); source("test-rescale-neighbors.R")

library(scater)    
exprs_sce <- mockSCE()
exprs_sce <- logNormCounts(exprs_sce)
exprs_sce <- runPCA(exprs_sce)

adt_sce <- mockSCE(ngenes=20) 
adt_sce <- logNormCounts(adt_sce)
altExp(exprs_sce, "ADT") <- adt_sce

test_that("rescaleByNeighbors works as expected", {
    combined <- rescaleByNeighbors(exprs_sce, dimreds="PCA", altexps="ADT")
    expect_identical(ncol(combined), nrow(adt_sce) + ncol(reducedDim(exprs_sce)))

    combined <- rescaleByNeighbors(exprs_sce, dimreds="PCA", altexps="ADT", combine=FALSE)
    expect_equal(1, median(BiocNeighbors::findKNN(combined[[1]], k=50, get.index=FALSE)$distance[,50]))
    expect_equal(1, median(BiocNeighbors::findKNN(combined[[2]], k=50, get.index=FALSE)$distance[,50]))

    combined2 <- rescaleByNeighbors(exprs_sce, dimreds="PCA", altexps="ADT", combine=FALSE)
    expect_true(sd(combined2[[1]]/reducedDim(exprs_sce)) < 1e-8)
    ratios <- combined2[[2]]/t(logcounts(adt_sce))
    ratios <- ratios[!is.nan(ratios)] # ignore zeros before and after.
    expect_true(sd(ratios) < 1e-8)
})

test_that("rescaleByNeighbors weighting is correct", {
    ref <- rescaleByNeighbors(exprs_sce, dimreds="PCA", altexps="ADT")

    # Robust to increases in the scope of one of the modes.
    alt <- exprs_sce
    logcounts(altExp(alt, "ADT")) <- 2*logcounts(altExp(alt, "ADT"))
    doubled <- rescaleByNeighbors(alt, dimreds="PCA", altexps="ADT")
    expect_equal(ref, doubled)

    alt <- exprs_sce
    altExp(alt, "ADT") <- rbind(altExp(alt, "ADT"), altExp(alt, "ADT"))
    doubled <- rescaleByNeighbors(alt, dimreds="PCA", altexps="ADT")
    expect_equal(as.matrix(dist(ref)), as.matrix(dist(doubled)))
})

test_that("rescaleByNeighbors works with mixed DelayedArrays and matrices", {
    ref <- rescaleByNeighbors(exprs_sce, dimreds="PCA", altexps="ADT")
    combined <- rescaleByNeighbors(list(reducedDim(exprs_sce), DelayedArray::DelayedArray(t(logcounts(adt_sce)))))
    expect_s4_class(combined, "DelayedArray")
    expect_identical(ref, as.matrix(combined))
})

test_that("rescaleByNeighbors custom weighting is correct", {
    out <- rescaleByNeighbors(exprs_sce, dimreds="PCA", altexps="ADT", weights=c(1,sqrt(2)))
    ref <- rescaleByNeighbors(exprs_sce, dimreds="PCA", altexps="ADT", extras=list(t(logcounts(altExp(exprs_sce, "ADT")))))
    expect_equal(as.matrix(dist(ref)), as.matrix(dist(out)))

    # Order of weights matches up with order of arguments for SCEs.
    test <- exprs_sce[1:10,]
    out <- rescaleByNeighbors(test, assays="logcounts", dimreds="PCA", altexps="ADT", weights=c(1,2,3), combine=FALSE)
    ref <- rescaleByNeighbors(list(t(logcounts(test)), reducedDim(test), t(logcounts(altExp(test, "ADT")))), weights=c(1,2,3), combine=FALSE)
    expect_identical(out, ref)
})

test_that("rescaleByNeighbors works correctly with SEs and extras", {
    test <- exprs_sce[1:10,]
    ref <- rescaleByNeighbors(test, assays="logcounts", dimreds="PCA", altexps="ADT")

    se <- as(test, "SummarizedExperiment")
    out <- rescaleByNeighbors(se, assays="logcounts", extras=list(reducedDim(test), t(logcounts(altExp(test, "ADT")))))

    expect_identical(unname(out), unname(ref))
})
