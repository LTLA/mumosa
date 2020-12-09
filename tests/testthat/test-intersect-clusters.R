# This tests the cluster intersection functionality.
# library(testthat); library(mumosa); source("test-intersect-clusters.R")

set.seed(10000)
mat1 <- matrix(rnorm(10000), ncol=20)
chosen <- 1:250
mat1[chosen,1] <- mat1[chosen,1] + 10
kout1 <- kmeans(mat1, 5)
clusters1 <- kout1$cluster

mat2 <- matrix(rnorm(10000), ncol=20)
chosen <- c(1:125, 251:375)
mat2[chosen,2] <- mat2[chosen,2] + 10
clusters2 <- kmeans(mat2, 5)$cluster

test_that("WCSS calculators work as expected", {
    expect_equal(
        mumosa:::.compute_wcss_solo(mat1),
        sum(t(t(mat1) - colMeans(mat1))^2)
    )

    expect_equal(
        unname(mumosa:::.compute_wcss(clusters1, mat1)),
        kout1$withinss
    )
})

test_that("intersectClusters works as expected", {
    clusters3 <- intersectClusters(list(clusters1, clusters2), list(mat1, mat2))
    expect_true(
        sum(mumosa:::.compute_wcss(clusters3, mat1)) < sum(mumosa:::.compute_wcss(clusters1, mat1))
    )
    expect_true(
        sum(mumosa:::.compute_wcss(clusters3, mat2)) < sum(mumosa:::.compute_wcss(clusters2, mat2))
    )
    
    tab <- table(clusters3, (mat1[,1] > 5) + 2*(mat2[,2] > 5))
    expect_identical(ncol(tab), 4L) # all four states present.
    expect_true(all(rowSums(tab > 0)==1)) # each cluster is only in one state.
})
