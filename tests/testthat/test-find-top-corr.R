# This tests the findTopCorrelations function.
# library(testthat); library(mumosa); source("test-find-top-corr.R")

set.seed(109109100)
test_that("findTopCorrelations works for the top self-correlations", {
    library(scuttle)
    sce <- mockSCE(ngenes=100)
    sce <- logNormCounts(sce)

    # Setting d to the max, so that the search is basically exact.
    df <- findTopCorrelations(sce, number=20, d=nrow(sce), BSPARAM=BiocSingular::ExactParam()) 
    df2 <- findTopCorrelations(sce, number=20, d=NA)
    expect_identical(df, df2)

    # Performing a reference computation.
    all.cor <- cor(t(logcounts(sce)), method="spearman")
    collected.best <- collected.worst <- collected.best.cor <- collected.worst.cor <- vector("list", nrow(all.cor))
    for (i in seq_len(nrow(all.cor))) {
        current <- setdiff(order(all.cor[i,], decreasing=TRUE), i)

        my.best <- head(current, 20)
        collected.best[[i]] <- colnames(all.cor)[my.best]
        collected.best.cor[[i]] <- all.cor[i,my.best]

        my.worst <- rev(tail(current, 20))
        collected.worst[[i]] <- colnames(all.cor)[my.worst]
        collected.worst.cor[[i]] <- all.cor[i,my.worst]
    }

    expect_identical(df$positive$gene2, unlist(collected.best))
    expect_equal(df$positive$rho, unname(unlist(collected.best.cor)))
    expect_identical(df$negative$gene2, unlist(collected.worst))
    expect_equal(df$negative$rho, unname(unlist(collected.worst.cor)))

    # Expect the same results from one setting of direction.
    dfp <- findTopCorrelations(sce, number=20, d=nrow(sce), direction="positive", BSPARAM=BiocSingular::ExactParam()) 
    expect_identical(dfp, df["positive"])

    dfn <- findTopCorrelations(sce, number=20, d=nrow(sce), direction="negative", BSPARAM=BiocSingular::ExactParam()) 
    expect_identical(dfn, df["negative"])
})

test_that("findTopCorrelations p-value is correctly computed", {
    mat <- matrix(runif(2000), ncol=100)
    df <- findTopCorrelations(mat, number=1, d=NA)

    up.p <- down.p <- numeric(nrow(df$positive))
    for (x in seq_len(nrow(df$positive))) {
        up.p[x] <- cor.test(mat[x,], mat[df$positive$gene2[x],], method="spearman", alternative="greater", exact=FALSE)$p.value
        down.p[x] <- cor.test(mat[x,], mat[df$negative$gene2[x],], method="spearman", alternative="less", exact=FALSE)$p.value
    }

    expect_equal(up.p, df$positive$p.value)
    expect_equal(down.p, df$negative$p.value)

    # This should not be true.
    expect_false(isTRUE(all.equal(p.adjust(df$p.value, method="BH"), df$FDR)))
})

test_that("findTopCorrelations correctly excludes self from negative correlations", {
    mat <- matrix(1:5, 10, 5, byrow=TRUE)
    df <- suppressWarnings(findTopCorrelations(mat, number=20, d=ncol(mat), BSPARAM=BiocSingular::ExactParam()) )

    expect_identical(nrow(df$positive), nrow(mat) * (nrow(mat)  - 1L))
    expect_identical(nrow(df$positive), nrow(df$negative))
    expect_true(all(df$negative$gene1!=df$negative$gene2))
})

set.seed(109109101)
test_that("findTopCorrelations approximations work well enough", {
    thingy <- matrix(runif(2000), ncol=100)
    mat <- rbind(thingy, thingy)

    df <- suppressWarnings(findTopCorrelations(mat, number=1, d=NA))
    set.seed(1000)
    df2 <- suppressWarnings(findTopCorrelations(mat, number=1, d=5))
    expect_identical(df$positive, df2$positive)

    mat <- rbind(thingy, -thingy)
    df <- suppressWarnings(findTopCorrelations(mat, number=1, d=NA))
    set.seed(1000)
    df2 <- suppressWarnings(findTopCorrelations(mat, number=1, d=5))
    expect_identical(df$negative, df2$negative)
})

set.seed(109109102)
test_that("findTopCorrelations works for the top cross-correlations", {
    sce1 <- mockSCE(ngenes=100)
    sce1 <- logNormCounts(sce1)

    sce2 <- mockSCE(ngenes=200)
    sce2 <- logNormCounts(sce2)

    df <- findTopCorrelations(sce1, y=sce2, number=20, d=ncol(sce1), BSPARAM=BiocSingular::ExactParam()) 
    df2 <- findTopCorrelations(sce1, y=sce2, number=20, d=NA)
    expect_identical(df, df2)

    # Performing a reference computation.
    all.cor <- cor(t(logcounts(sce1)), t(logcounts(sce2)), method="spearman")
    collected.best <- collected.worst <- collected.best.cor <- collected.worst.cor <- vector("list", nrow(all.cor))
    for (i in seq_len(nrow(all.cor))) {
        current <- order(all.cor[i,], decreasing=TRUE)

        my.best <- head(current, 20)
        collected.best[[i]] <- colnames(all.cor)[my.best]
        collected.best.cor[[i]] <- all.cor[i,my.best]

        my.worst <- rev(tail(current, 20))
        collected.worst[[i]] <- colnames(all.cor)[my.worst]
        collected.worst.cor[[i]] <- all.cor[i,my.worst]
    }

    expect_identical(df$positive$gene2, unlist(collected.best))
    expect_equal(df$positive$rho, unname(unlist(collected.best.cor)))
    expect_identical(df$negative$gene2, unlist(collected.worst))
    expect_equal(df$negative$rho, unname(unlist(collected.worst.cor)))

    # Expect the same results from one setting of direction.
    dfp <- findTopCorrelations(sce1, y=sce2, number=20, d=ncol(sce1), direction="positive", BSPARAM=BiocSingular::ExactParam()) 
    expect_identical(dfp, df["positive"])

    dfn <- findTopCorrelations(sce1, y=sce2, number=20, d=ncol(sce2), direction="negative", BSPARAM=BiocSingular::ExactParam()) 
    expect_identical(dfn, df["negative"])
})

set.seed(109109102)
test_that("findTopCorrelations works with blocked self-correlations", {
    library(scuttle)
    sce <- mockSCE(ngenes=100)
    sce <- logNormCounts(sce)

    library(mumosa)
    ref <- findTopCorrelations(sce, number=20, d=10, BSPARAM=BiocSingular::ExactParam())
    out <- findTopCorrelations(sce, number=20, d=10, block=rep(1, ncol(sce)), BSPARAM=BiocSingular::ExactParam()) 
    expect_identical(ref,out)

    # Equiweighting preserves the results when blocks are unbalanced.
    block <- rep(1:2, each=ncol(sce)/2)
    ref <- findTopCorrelations(sce, number=20, d=10, block=block, BSPARAM=BiocSingular::ExactParam())

    expanded <- c(seq_len(ncol(sce)), seq_len(ncol(sce)/2))
    out <- findTopCorrelations(sce[,expanded], number=20, d=10, block=block[expanded], BSPARAM=BiocSingular::ExactParam()) 

    ref$positive <- ref$positive[,1:3]
    ref$positive <- ref$positive[do.call(order, ref$positive),]
    out$positive <- out$positive[,1:3]
    out$positive <- out$positive[do.call(order, out$positive),]
    ref$negative <- ref$negative[,1:3]
    ref$negative <- ref$negative[do.call(order, ref$negative),]
    out$negative <- out$negative[,1:3]
    out$negative <- out$negative[do.call(order, out$negative),]
    expect_equal(ref, out)

    # Same results without equiweighting, for balanced blocks.
    weight <- findTopCorrelations(sce, number=20, d=NA, block=block, BSPARAM=BiocSingular::ExactParam())
    noweight <- findTopCorrelations(sce, number=20, d=NA, block=block, equiweight=FALSE, BSPARAM=BiocSingular::ExactParam())
    expect_equal(weight, noweight)

    # Not true of unweighted and unbalanced blocks.    
    block0 <- rep(1:2, c(20, ncol(sce)-20))
    weight <- findTopCorrelations(sce, number=20, d=NA, block=block0, BSPARAM=BiocSingular::ExactParam())
    noweight <- findTopCorrelations(sce, number=20, d=NA, block=block0, equiweight=FALSE, BSPARAM=BiocSingular::ExactParam())
    expect_false(isTRUE(all.equal(weight, noweight)))
})



