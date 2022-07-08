context("Test inference under no bias")

test_that("Selected bw is infinite", {
    d <- FRDData(cbind(logcn=log(rcp[1:5000, 6]), rcp[1:5000, c(3, 2)]),
                 cutoff=0)

    ## Expect using all data
    r0 <- NPRHonest.fit(d, M=c(0, 0), kern="triangular", opt.criterion="OCI",
                        T0=0)$coefficients
    r1 <- NPRHonest.fit(d, M=c(0, 0), kern="uniform", opt.criterion="FLCI",
                        T0=0)$coefficients
    r2 <- NPRreg.fit(d, Inf, "uniform")
    expect_equal(c(r1$std.error, r1$estimate),
                 c(unname(r2$se["nn"]), r2$estimate))
    expect_identical(r1$maximum.bias, 0)
    ## For triangular kernel, maximum bw is given by end of support, so we
    ## expect to find minimum at boundary
    expect_equal(r1$bandwidth, max(abs(c(d$Xp, d$Xm))))

    d <- RDData(lee08, cutoff=0)
    r1 <- NPRHonest.fit(d, M=0, kern="uniform",
                        opt.criterion="MSE")$coefficients
    r2 <- NPRreg.fit(d, Inf, "uniform")
    expect_equal(c(r1$std.error, r1$estimate),
                 c(unname(r2$se["nn"]), r2$estimate))
    expect_identical(r1$maximum.bias, 0)
})

context("Test FRD")

test_that("FRD data example check", {
    d <- FRDData(cbind(logcn=log(rcp[, 6]), rcp[, c(3, 2)]), cutoff=0)
    M <- NPR_MROT.fit(d)
    r1 <- NPRHonest.fit(d, M, kern="triangular",
                        opt.criterion="MSE", T0=0)$coefficients
    r2 <- NPRHonest.fit(d, M, kern="triangular", opt.criterion="MSE",
                        T0=r1$estimate)$coefficients
    ## With positive T0, expect greater effective M, and thus smaller bandwidth
    expect_lt(r2$bandwidth, r1$bandwidth)
    ## On travis, only the first 6 digits match, not sure why
    expect_equal(round(r2$estimate, 6), -0.188134)
})

test_that("FRD with almost perfect first stage", {
    d <- FRDData(data.frame(y=lee08$voteshare,
                            d=lee08$margin>0, x=lee08$margin), cutoff=0)
    M <- NPR_MROT.fit(d)
    r0 <- NPRHonest.fit(d, M, kern="triangular", opt.criterion="FLCI")
    r1 <- NPRHonest.fit(d, M, kern="triangular", opt.criterion="FLCI",
                        T0=r0$coefficients$estimate)$coefficients
    r2 <- NPRHonest.fit(RDData(lee08, cutoff=0), unname(M[1]),
                        kern="triangular", opt.criterion="FLCI")$coefficients
    expect_lt(max(abs(r2[2:6]-r2[2:6])), 3e-7)

    df <- data.frame(y=lee08$voteshare,
                     d=lee08$margin+rnorm(n=length(lee08$margin), sd=0.1)>0,
                     x=lee08$margin)
    d <- FRDData(df, cutoff=0)
    M <- NPR_MROT.fit(d)
    r0 <- NPRHonest.fit(d, M, kern="triangular", opt.criterion="MSE")
    r1 <- NPRHonest.fit(d, M, kern="triangular", opt.criterion="MSE",
                        T0=r0$coefficients$estimate)
    r2 <- NPRHonest.fit(RDData(lee08, cutoff=0), unname(M[1]),
                        kern="triangular", opt.criterion="MSE")$coefficients
    expect_lt(abs(r1$coefficients$estimate*r1$fs-r2$estimate), 1e-2)
    expect_lt(abs(r1$coefficients$bandwidth-r2$bandwidth), 0.1)
    expect_lt(abs(r0$coefficients$conf.low-r1$coefficients$conf.low), 0.01)
})

test_that("FRD interface", {
    d <- FRDData(cbind(logf=log(rcp[1:10000, 6]), rcp[1:10000, c(3, 2)]),
                 cutoff=0)
    M <- NPR_MROT.fit(d)
    rcp1 <- rcp[1:10000, ]
    r1 <- NPRHonest.fit(d, M, kern="triangular", opt.criterion="OCI", T0=0)
    p1 <- RDHonest(log(cn)~retired | elig_year, data=rcp1, cutoff=0, M=M,
                    kern="triangular", opt.criterion="OCI", T0=0)
    expect_equal(r1$coefficients$estimate, p1$coefficients$estimate)

    r2 <- NPRHonest.fit(d, M, kern="triangular", opt.criterion="OCI",
                        T0=r1$coefficients$estimate)
    p2 <- RDHonest(log(cn)~retired | elig_year, data=rcp1, cutoff=0, M=M,
                    kern="triangular", opt.criterion="OCI",
                    T0=p1$coefficients$estimate)
    expect_equal(r2$coefficients$estimate, p2$coefficients$estimate)
    ## codecov.io check
    expect_equal(capture.output(print(r2))[1:7],
                 capture.output(print(p2))[7:13])

    r3 <- NPRHonest.fit(d, M, kern="triangular", h=7,
                        T0=r1$coefficients$estimate)
    p3 <- RDHonest(log(cn)~retired | elig_year, data=rcp1, cutoff=0, M=M,
                    kern="triangular", opt.criterion="OCI",
                    T0=p1$coefficients$estimate, h=7)
    expect_equal(r3$coefficients$estimate, p3$coefficients$estimate)

    r4 <- NPROptBW.fit(d, M, kern="triangular", opt.criterion="OCI",
                        T0=r1$coefficients$estimate)
    p4 <- RDHonest(log(cn)~retired | elig_year, data=rcp1, cutoff=0, M=M,
                    kern="triangular", opt.criterion="OCI",
                    T0=p1$coefficients$estimate)
    expect_identical(r4$h, p4$coefficients$bandwidth)
})
