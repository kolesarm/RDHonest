test_that("Test weighting using cghs", {
    ## Make 10 groups, estimate within-group variance
    set.seed(42)
    ## Make cells by group and year
    cghs$cell <- 0.1*sample((1:10), size=NROW(cghs), replace=TRUE) +
        cghs$yearat14

    dd <- data.frame()
    for (j in sort(unique(cghs$cell))) {
        ix <- cghs$cell==j
        xj <- cghs$yearat14[which.max(ix)]
        df <- data.frame(y=mean(log(cghs$earnings[ix])), x=xj-1947,
                         weights=sum(ix),
                         sigma2=var(log(cghs$earnings[cghs$yearat14==xj])) /
                             sum(ix))
        dd <- rbind(dd, df)
    }
    s0 <- RDHonest(log(earnings)~yearat14, cutoff=1947, h=5,
                   data=cghs, M=1)$coefficients
    s1 <- RDHonest(y~x, h=5, data=dd, M=1, weights=weights,
                   sigmaY2=sigma2, se.method="supplied.var")$coefficients
    rownames(s0) <- rownames(s1)
    expect_equal(s0, s1)

    ## Don't supply variance
    s2 <- RDHonest(y~x, h=5, data=dd, M=1, weights=weights)$coefficients
    expect_lt(max(abs(s0[c(2, 4, 10, 11)] - s2[c(2, 4, 10, 11)])), 1e-8)
    ## SE should be close enough
    expect_lt(max(abs(s2[c(3, 5, 6)]-s0[c(3, 5, 6)])), 3e-3)
    ## Estimate M
    expect_message(s0a <- RDHonest(log(earnings)~yearat14, cutoff=1947, h=5,
                                   data=cghs, kern="uniform")$coefficients)
    expect_message(s2a <- RDHonest(y~x, h=5, data=dd, weights=weights,
                                   kern="uniform")$coefficients)
    expect_lt(max(abs(s0a[c(2, 4, 9:11)]-s2a[c(2, 4, 9:11)])), 1e-8)
    ## SE should be close enough
    expect_lt(max(abs(s2a[c(3, 5, 6)] - s0a[c(3, 5, 6)])), 3e-3)

    ## LPP Honest
    t0 <- RDHonest(log(earnings)~yearat14, cutoff=1947, h=5,
                   data=cghs[cghs$yearat14>=1947, ], M=1,
                   point.inference=TRUE)$coefficients
    t1 <- RDHonest(y~x, cutoff=0, h=5, data=dd[dd$x>=0, ], M=1,
                   weights=weights, point.inference=TRUE,
                   sigmaY2=sigma2, se.method="supplied.var")$coefficients
    t2 <- RDHonest(y~x, cutoff=0, h=5, data=dd[dd$x>=0, ], M=1,
                   weights=weights, point.inference=TRUE)$coefficients
    expect_equal(t0, t1)
    expect_lt(max(abs(t0[c(2, 4, 9:11)] - t2[c(2, 4, 9:11)])), 1e-8)
    expect_lt(max(abs(t2[c(3, 5, 6)]-t0[c(3, 5, 6)])), 1.1e-3)

    ## Collapse data fully
    expect_message(s0 <- RDHonest(log(earnings)~yearat14,
                                  cutoff=1947, data=cghs)$coefficients)
    dd <- data.frame()
    for (j in unique(cghs$yearat14)) {
        ix <- cghs$yearat14==j
        df <- data.frame(y=mean(log(cghs$earnings[ix])), x=j,
                         weights=sum(ix),
                         sigma2=var(log(cghs$earnings[ix]))/sum(ix))
        dd <- rbind(dd, df)
    }
    expect_message(s1 <- RDHonest(y~x, cutoff=1947, data=dd, weights=weights,
                                  sigmaY2=sigma2, se.method="supplied.var",
                                  h=s0$bandwidth))
    expect_equal(rownames(s0), "I(yearat14>0)")
    rownames(s0) <- rownames(s1$coefficients)
    expect_equal(s1$coefficients, s0)
    ## If we use supplied.var for variance estimation, should match approx
    expect_message(s2 <- RDHonest(y~x, cutoff=1947, data=dd,
                                  weights=weights, sigmaY2=sigma2,
                                  se.method="supplied.var")$coefficients)
    expect_lt(max(abs(s2[2:9]-s0[2:9])), 4e-3)

    ## Test that sigmaNN works when J large
    expect_message(s3 <- RDHonest(y~x, cutoff=1947, data=dd, weights=weights,
                                  J=10, h=s0$bandwidth)$coefficients)
    expect_message(s4 <- RDHonest(y~x, cutoff=1947, data=dd, weights=weights,
                                  J=20, h=s0$bandwidth)$coefficients)
    expect_equal(s3, s4)
})
## Fuzzy RD
test_that("Test weighting using rcp ", {
    expect_message(r0 <- RDHonest(log(cn)|retired~elig_year,
                                  data=rcp, T0=0)$coefficients)
    dd <- data.frame(x=sort(unique(rcp$elig_year)), y=NA, d=NA, weights=NA,
                     sig11=NA, sig12=NA, sig21=NA, sig22=NA)
    for (j in seq_len(NROW(dd))) {
        ix <- rcp$elig_year==dd$x[j]
        Y <- cbind(log(rcp$cn[ix]), rcp$retired[ix])
        dd[j, -1] <- c(colMeans(Y), sum(ix), as.vector(var(Y))/sum(ix))
    }
    expect_message(r1 <- RDHonest(y|d~x, data=dd, weights=weights,
                                  sigmaY2=sig11, T0=0, sigmaYD=sig21,
                                  sigmaD2=sig22, h=r0$bandwidth,
                                  se.method="supplied.var")$coefficients)
    expect_equal(rownames(r0), "retired")
    rownames(r1) <- rownames(r0)
    expect_equal(r0, r1)
    expect_message(r2 <- RDHonest(y|d~x, data=dd, weights=weights, T0=0,
                                  sigmaY2=sig11, sigmaYD=sig21, sigmaD2=sig22,
                                  se.method="supplied.var")$coefficients)
    expect_lt(max(abs(r2[2:8]-r0[2:8])), 2e-3)
})
