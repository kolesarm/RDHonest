context("Test inference under no bias")

test_that("Selected bw is infinite", {
    d <- FRDData(cbind(logcn=log(rcp[1:5000, 6 ]), rcp[1:5000, c(3, 2)]), cutoff=0)

    ## Expect using all data
    r0 <- NPRHonest.fit(d, M=c(0, 0), kern="triangular",
                        opt.criterion="OCI", T0=0)
    r1 <- NPRHonest.fit(d, M=c(0, 0), kern="uniform",
                        opt.criterion="FLCI", T0=0)
    r2 <- NPRreg.fit(d, h=Inf, kern="uniform")
    expect_equal(c(r1$sd, r1$estimate), c(r2$se["nn"], r2$estimate))
    expect_identical(r1$maxbias, 0)
    ## For triangular kernel, maximum bw is given by end of support, so we
    ## expect to find minimum at boundary
    expect_equal(r1$hp, max(abs(c(d$Xp, d$Xm))))

    d <- RDData(lee08, cutoff=0)
    r1 <- NPRHonest.fit(d, M=0, kern="uniform", opt.criterion="MSE")
    r2 <- NPRreg.fit(d, h=Inf, kern="uniform")
    expect_equal(c(r1$sd, r1$estimate), c(r2$se["nn"], r2$estimate))
    expect_identical(r1$maxbias, 0)
})

context("Test FRD")

test_that("FRD data example check", {

    d <- FRDData(cbind(logcn=log(rcp[, 6 ]), rcp[, c(3, 2)]), cutoff=0)
    M <- NPR_MROT.fit(d)
    r1 <- NPRHonest.fit(d, M, kern="triangular", opt.criterion="MSE", T0=0)
    r2 <- NPRHonest.fit(d, M, kern="triangular", opt.criterion="MSE",
                        T0=r1$estimate)
    ## With positive T0, expect greater effective M, and thus smaller bandwidth
    expect_lt(r2$hp, r1$hp)
    expect_equal(r2$estimate, -0.18813375)
})

test_that("FRD with almost perfect first stage", {
    d <- FRDData(data.frame(y=lee08$voteshare,
                            d=lee08$margin>0, x=lee08$margin), cutoff=0)
    M <- NPR_MROT.fit(d)
    r0 <- NPRHonest.fit(d, M, kern="triangular", opt.criterion="FLCI")
    r1 <- NPRHonest.fit(d, M, kern="triangular", opt.criterion="FLCI",
                        T0=r0$estimate)
    r2 <- NPRHonest.fit(RDData(lee08, cutoff=0), unname(M[1]),
                        kern="triangular", opt.criterion="FLCI")
    expect_equal(c(r1$estimate, r1$sd, r1$hl),
                 c(r2$estimate, r2$sd, r2$hl))

    df <- data.frame(y=lee08$voteshare,
                     d=lee08$margin+rnorm(n=length(lee08$margin), sd=0.1)>0,
                     x=lee08$margin)
    d <- FRDData(df, cutoff=0)
    M <- NPR_MROT.fit(d)
    r0 <- NPRHonest.fit(d, M, kern="triangular", opt.criterion="MSE")
    r1 <- NPRHonest.fit(d, M, kern="triangular", opt.criterion="MSE",
                        T0=r0$estimate)
    r2 <- NPRHonest.fit(RDData(lee08, cutoff=0), unname(M[1]),
                        kern="triangular", opt.criterion="MSE")
    expect_lt(abs(r1$estimate*r1$fs-r2$estimate), 1e-2)
    expect_lt(abs(r1$hp-r2$hp), 0.1)
    expect_lt(abs(r0$hl-r1$hl), 0.01)
})

test_that("FRD interface", {
    d <- FRDData(cbind(logf=log(rcp[1:10000, 6]), rcp[1:10000, c(3, 2)]), cutoff=0)
    M <- NPR_MROT.fit(d)
    rcp1 <- rcp[1:10000, ]
    r1 <- NPRHonest.fit(d, M, kern="triangular", opt.criterion="OCI", T0=0)
    p1 <- FRDHonest(log(cn)~retired | elig_year, data=rcp1, cutoff=0, M=M,
                    kern="triangular", opt.criterion="OCI", T0=0)
    expect_equal(r1$estimate, p1$estimate)

    r2 <- NPRHonest.fit(d, M, kern="triangular", opt.criterion="OCI",
                        T0=r1$estimate)
    p2 <- FRDHonest(log(cn)~retired | elig_year, data=rcp1, cutoff=0, M=M,
                    kern="triangular", opt.criterion="OCI", T0=p1$estimate)
    expect_equal(r2$estimate, p2$estimate)
    expect_equal(capture.output(print(r2)), capture.output(print(p2))[-(1:6)])

    r3 <- NPRHonest.fit(d, M, kern="triangular", h=7,
                        T0=r1$estimate)
    p3 <- FRDHonest(log(cn)~retired | elig_year, data=rcp1, cutoff=0, M=M,
                    kern="triangular", opt.criterion="OCI", T0=p1$estimate, h=7)
    expect_equal(r3$estimate, p3$estimate)

    r4 <- NPROptBW.fit(d, M, kern="triangular", opt.criterion="OCI",
                        T0=r1$estimate)
    p4 <- FRDOptBW(log(cn)~retired | elig_year, data=rcp1, cutoff=0, M=M,
                    kern="triangular", opt.criterion="OCI", T0=p1$estimate)
    expect_equal(capture.output(print(r4)), capture.output(print(p4))[7])
})
