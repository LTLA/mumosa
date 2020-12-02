# This tests the findTopCorrelations function.
# library(testthat); library(mumosa); source("setup.R"); source("test-find-top-corr.R")

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

    expect_identical(df$positive$feature2, unlist(collected.best))
    expect_equal(df$positive$rho, unname(unlist(collected.best.cor)))
    expect_identical(df$negative$feature2, unlist(collected.worst))
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
        up.p[x] <- cor.test(mat[x,], mat[df$positive$feature2[x],], method="spearman", alternative="greater", exact=FALSE)$p.value
        down.p[x] <- cor.test(mat[x,], mat[df$negative$feature2[x],], method="spearman", alternative="less", exact=FALSE)$p.value
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
    expect_true(all(df$negative$feature1!=df$negative$feature2))
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

    expect_identical(df$positive$feature2, unlist(collected.best))
    expect_equal(df$positive$rho, unname(unlist(collected.best.cor)))
    expect_identical(df$negative$feature2, unlist(collected.worst))
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

set.seed(109109103)
test_that("findTopCorrelations works with blocked cross-correlations", {
    sce1 <- mockSCE(ngenes=100)
    sce1 <- logNormCounts(sce1)

    sce2 <- mockSCE(ngenes=200)
    sce2 <- logNormCounts(sce2)

    ref <- findTopCorrelations(sce1, y=sce2, number=20, d=10, BSPARAM=BiocSingular::ExactParam())
    out <- findTopCorrelations(sce1, y=sce2, number=20, d=10, block=rep(1, ncol(sce1)), BSPARAM=BiocSingular::ExactParam()) 
    expect_identical(ref,out)

    # Equiweighting preserves the results when blocks are unbalanced.
    block <- rep(1:2, each=ncol(sce1)/2)
    ref <- findTopCorrelations(sce1, y=sce2, number=20, d=10, block=block, BSPARAM=BiocSingular::ExactParam())

    expanded <- c(seq_len(ncol(sce1)), seq_len(ncol(sce1)/2))
    out <- findTopCorrelations(sce1[,expanded], y=sce2[,expanded], number=20, d=10, 
        block=block[expanded], BSPARAM=BiocSingular::ExactParam()) 

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
    weight <- findTopCorrelations(sce1, y=sce2, number=20, d=NA, block=block, BSPARAM=BiocSingular::ExactParam())
    noweight <- findTopCorrelations(sce1, y=sce2, number=20, d=NA, block=block, equiweight=FALSE, BSPARAM=BiocSingular::ExactParam())
    expect_equal(weight, noweight)

    # Not true of unweighted and unbalanced blocks.    
    block0 <- rep(1:2, c(20, ncol(sce1)-20))
    weight <- findTopCorrelations(sce1, y=sce2, number=20, d=NA, block=block0, BSPARAM=BiocSingular::ExactParam())
    noweight <- findTopCorrelations(sce1, y=sce2, number=20, d=NA, block=block0, equiweight=FALSE, BSPARAM=BiocSingular::ExactParam())
    expect_false(isTRUE(all.equal(weight, noweight)))
})

set.seed(109109104)
test_that("findTopCorrelations works with unblocked zero-variance genes", {
    sce1 <- mockSCE(ngenes=100)
    vals1 <- counts(sce1)
    vals1[1:10,] <- 0:9

    ref <- findTopCorrelations(vals1, number=20, d=10, BSPARAM=BiocSingular::ExactParam())
    out <- findTopCorrelations(vals1[-(1:10),], number=20, d=10, BSPARAM=BiocSingular::ExactParam())
    expect_identical(ref, out)

    sce2 <- mockSCE(ngenes=100)
    vals2 <- counts(sce2)
    vals2[1:10,] <- 0:9

    ref <- findTopCorrelations(vals1, y=vals2, number=20, d=10, BSPARAM=BiocSingular::ExactParam())
    out <- findTopCorrelations(vals1[-(1:10),], y=vals2[-(1:10),], number=20, d=10, BSPARAM=BiocSingular::ExactParam())
    expect_identical(ref, out)
})

set.seed(109109105)
test_that("findTopCorrelations (self) works with blocked zero-variance genes", {
    sce1 <- mockSCE(ngenes=100)
    vals1 <- counts(sce1)

    block <- sample(2, ncol(sce1), replace=TRUE)
    vals1[1:5,block==1] <- 0
    vals1[6:10,block==2] <- 0

    out <- suppressWarnings(findTopCorrelations(vals1[1:10,], number=20, d=10, block=block, BSPARAM=BiocSingular::ExactParam()))
    ref1 <- suppressWarnings(findTopCorrelations(vals1[1:5,block!=1], number=20, d=10, BSPARAM=BiocSingular::ExactParam()))
    ref2 <- suppressWarnings(findTopCorrelations(vals1[6:10,block!=2], number=20, d=10, BSPARAM=BiocSingular::ExactParam()))
    expect_identical(out$positive[,1:4], rbind(ref1$positive[,1:4], ref2$positive[,1:4]))
    expect_identical(out$negative[,1:4], rbind(ref1$negative[,1:4], ref2$negative[,1:4]))

    others <- 11:nrow(vals1)
    refo <- findTopCorrelations(vals1[others,], number=20, d=10, block=block, BSPARAM=BiocSingular::ExactParam())
    ref1o <- findTopCorrelations(vals1[1:5,block!=1], y=vals1[others,block!=1], number=20, d=10, BSPARAM=BiocSingular::ExactParam())
    refo1 <- suppressWarnings(findTopCorrelations(vals1[others,block!=1], y=vals1[1:5,block!=1], number=20, d=10, BSPARAM=BiocSingular::ExactParam()))
    ref2o <- findTopCorrelations(vals1[6:10,block!=2], y=vals1[others,block!=2], number=20, d=10, BSPARAM=BiocSingular::ExactParam())
    refo2 <- suppressWarnings(findTopCorrelations(vals1[others,block!=2], y=vals1[6:10,block!=2], number=20, d=10, BSPARAM=BiocSingular::ExactParam()))

    refp <- rbind(ref1$positive, ref2$positive, refo$positive, ref1o$positive, refo1$positive, ref2o$positive, refo2$positive)
    refp <- refp[order(refp$feature1, refp$p.value),]
    refp <- unlist(heads(split(refp, refp$feature1), 20), use.names=FALSE)
    refn <- rbind(ref1$negative, ref2$negative, refo$negative, ref1o$negative, refo1$negative, ref2o$negative, refo2$negative)
    refn <- refn[order(refn$feature1, refn$p.value),]
    refn <- unlist(heads(split(refn, refn$feature1), 20), use.names=FALSE)

    out <- suppressWarnings(findTopCorrelations(vals1, number=20, d=10, block=block, BSPARAM=BiocSingular::ExactParam()))
    expect_identical(out$positive[,1:4], refp[,1:4])
    expect_identical(out$negative[,1:4], refn[,1:4])
})

set.seed(109109106)
test_that("findTopCorrelations (self) works with blocked zero-variance genes", {
    library(scuttle)
    sce1 <- mockSCE(ngenes=100)
    vals1 <- counts(sce1)
    sce2 <- mockSCE(ngenes=200)
    vals2 <- counts(sce2)

    block <- sample(2, ncol(sce1), replace=TRUE)
    vals1[1:5,block==1] <- 0
    vals2[6:10,block==2] <- 0

    out <- suppressWarnings(findTopCorrelations(vals1[1:5,], y=vals2[6:10,], number=20, d=10, block=block, BSPARAM=BiocSingular::ExactParam()))
    expect_identical(nrow(out$positive), 0L)
    expect_identical(nrow(out$negative), 0L)

    others1 <- seq_len(nrow(vals1))[-(1:5)]
    others2 <- seq_len(nrow(vals2))[-(6:10)]
    refo <- findTopCorrelations(vals1[others1,], y=vals2[others2,], number=20, d=10, block=block, BSPARAM=BiocSingular::ExactParam())
    ref1 <- findTopCorrelations(vals1[1:5,block!=1], y=vals2[others2,block!=1], number=20, d=10, BSPARAM=BiocSingular::ExactParam())
    ref2 <- suppressWarnings(findTopCorrelations(vals1[others1,block!=2], y=vals2[6:10,block!=2], number=20, d=10, BSPARAM=BiocSingular::ExactParam()))

    refp <- rbind(ref1$positive, ref2$positive, refo$positive)
    refp <- refp[order(refp$feature1, refp$p.value),]
    refp <- unlist(heads(split(refp, refp$feature1), 20), use.names=FALSE)
    refn <- rbind(ref1$negative, ref2$negative, refo$negative)
    refn <- refn[order(refn$feature1, refn$p.value),]
    refn <- unlist(heads(split(refn, refn$feature1), 20), use.names=FALSE)

    out <- suppressWarnings(findTopCorrelations(vals1, y=vals2, number=20, d=10, block=block, BSPARAM=BiocSingular::ExactParam()))
    expect_identical(out$positive[,1:4], refp[,1:4])
    expect_identical(out$negative[,1:4], refn[,1:4])
})
