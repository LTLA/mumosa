# This tests the multi-UMAP code.
# library(testthat); library(mumosa); source("test-multi-umap.R")

stuff <- matrix(rnorm(10000), ncol=50)
things <- list(stuff, stuff[,1:5], stuff[,1:20])

test_that("metrics assembly works as expected", {
    metrics <- mumosa:::.compute_multi_modal_metrics(things)
    expect_equal(names(metrics), rep("euclidean", length(things)))
    expect_identical(unname(lengths(metrics)), vapply(things, ncol, 0L))
    expect_identical(unname(unlist(metrics)), seq_len(sum(vapply(things, ncol, 0L))))
})

test_that("multi-modal UMAP works as expected", {
    output <- calculateMultiUMAP(things)
    expect_identical(nrow(output), nrow(stuff))
    expect_identical(ncol(output), 2L)

    set.seed(9999)
    output <- calculateMultiUMAP(things, n_components=10)
    expect_identical(nrow(output), nrow(stuff))
    expect_identical(ncol(output), 10L)

    # Same result for SCEs.
    sce <- SingleCellExperiment(list(X=t(stuff)), reducedDims=list(Y=stuff[,1:5]), altExps=list(Z=SummarizedExperiment(t(stuff[,1:20]))))

    set.seed(9999)
    output2 <- runMultiUMAP(sce, assays=1, dimreds=1, altexps=1, altexp.assay=1, n_components=10)
    expect_identical(output, reducedDim(output2, "MultiUMAP"))
})

test_that("multi-modal UMAP works with mixed DelayedArrays and matrices", {
    set.seed(9999)
    output1 <- calculateMultiUMAP(things)

    set.seed(9999)
    output2 <- calculateMultiUMAP(lapply(things, DelayedArray::DelayedArray))
    expect_identical(output1, output2)

    # Same result for SCEs.
    sce <- SingleCellExperiment(list(X=DelayedArray::DelayedArray(t(stuff))),
                                reducedDims=list(Y=stuff[,1:5]), altExps=list(Z=SummarizedExperiment(t(stuff[,1:20]))))

    set.seed(9999)
    output3 <- calculateMultiUMAP(things, n_components=10)

    set.seed(9999)
    output4 <- runMultiUMAP(sce, assays=1, dimreds=1, altexps=1, altexp.assay=1, n_components=10)
    expect_identical(output3, reducedDim(output4, "MultiUMAP"))
})
