test_that("Selected bw is infinite", {
    r0 <- RDHonest(log(cn)~retired|elig_year, data=rcp[1:5000, ], M=c(0, 0),
                   opt.criterion="OCI")
    ## For triangular kernel, maximum bw is given by end of support, so we
    ## expect to find minimum at boundary
    expect_gt(r0$coefficients$bandwidth, max(abs(r0$d$X))-2e-6)
    r0$d$sigma2 <- NULL
    r1 <- NPRHonest(r0$d, M=c(0, 0), kern="uniform", opt.criterion="FLCI",
                    T0=0)$coefficients
    r2 <- NPReg(r0$d, Inf, "uniform")
    expect_identical(r1$maximum.bias, 0)

    expect_equal(c(r1$std.error, r1$estimate),
                 unname(c(r2$se, r2$estimate)))
    expect_equal(r1$bandwidth, max(abs(r0$d$X)))

    r1 <- RDHonest(voteshare~margin, data=lee08, M=0, kern="uniform",
                   opt.criterion="MSE")
    r2 <- NPReg(r1$d, Inf, "uniform")
    expect_equal(c(r1$coefficients$std.error, r1$coefficients$estimate),
                 unname(c(r2$se, r2$estimate)))
    expect_identical(r1$coefficients$maximum.bias, 0)
})

test_that("FRD data example check", {
    expect_message(r1 <- RDHonest(log(cn)~retired|elig_year,
                                  data=rcp[1:5000, ]))
    r1$d$sigma2 <- NULL
    M0 <- c(r1$coefficients$M.rf, r1$coefficients$M.fs)
    r2 <- NPRHonest(r1$d, M=M0, kern="triangular", opt.criterion="MSE",
                    T0=r1$coefficients$estimate)$coefficients
    r1a <- NPRHonest(r1$d, M=M0, kern="triangular", opt.criterion="MSE", T0=0,
                     alpha=r1$coefficients$p.value)$coefficients

    expect_equal(r1a$conf.high, 0)
    ## With positive T0, expect greater effective M, and thus smaller bandwidth
    expect_lt(r2$bandwidth, r1$coefficients$bandwidth)
    ## On travis, only the first 6 digits match, not sure why
    expect_equal(round(r2$estimate, 6),  -0.244582)
})

test_that("FRD with almost perfect first stage", {
    expect_message(r0 <- RDHonest(voteshare~I(margin>0)|margin, data=lee08,
                                  kern="triangular", opt.criterion="FLCI"))
    r0$d$sigma2 <- NULL
    M0 <- c(r0$coefficients$M.rf, r0$coefficients$M.fs)
    r1 <- NPRHonest(r0$d, M0, kern="triangular", opt.criterion="FLCI",
                    T0=r0$coefficients$estimate)$coefficients
    r2 <- RDHonest(voteshare~margin, M=M0[1], data=lee08,
                   kern="triangular", opt.criterion="FLCI")$coefficients
    expect_lt(max(abs(r1[2:8]-r2[2:8])), 3e-7)

    set.seed(42)
    df <- data.frame(y=lee08$voteshare,
                     d=lee08$margin+rnorm(n=length(lee08$margin), sd=0.1)>0,
                     x=lee08$margin)
    expect_message(r0 <- RDHonest(y~d|x, data=df,
                                  kern="triangular", opt.criterion="MSE"))
    r0$d$sigma2 <- NULL
    M0 <- c(r0$coefficients$M.rf, r0$coefficients$M.fs)
    r1 <- NPRHonest(r0$d, M0, kern="triangular", opt.criterion="MSE",
                    T0=r0$coefficients$estimate)$coefficients
    r2 <- RDHonest(voteshare~margin, data=lee08, M=M0[1], kern="triangular",
                   opt.criterion="MSE")$coefficients
    expect_lt(abs(r1$estimate*r1$first.stage-r2$estimate), 2e-3)
    expect_lt(abs(r1$bandwidth-r2$bandwidth), 2e-2)
    expect_lt(abs(r0$coefficients$conf.low-r1$conf.low), 1e-3)
})

test_that("FRD interface", {
    rcp1 <- rcp[1:10000, ]
    expect_message(p1.1 <- RDHonest(log(cn)~retired | elig_year,
                                    data=rcp1, cutoff=0,
                                    kern="triangular", T0=0, h=6))
    expect_lt(abs(p1.1$coefficients$conf.high-0.150118293818), 1e-11)
    expect_lt(abs(p1.1$coefficients$maximum.bias-0.0980630668382), 1e-11)
    expect_message(p1.2 <- RDHonest(log(cn)~retired | elig_year,
                                    data=rcp1, cutoff=0))
    expect_lt(abs(p1.2$coefficients$M.fs-0.00626732839127), 1e-11)
    expect_lt(abs(p1.2$coefficients$M.rf-0.00436022188088), 1e-11)
    expect_lt(abs(p1.2$coefficients$conf.high- 0.17195578263), 1e-6)
    expect_lt(abs(p1.2$coefficients$estimate+0.172378084757), 1e-6)

    p2 <- RDHonest(log(cn)~retired | elig_year, data=rcp1, cutoff=0,
                   M=as.numeric(p1.1$coefficients[16:17]),
                   kern="triangular", opt.criterion="OCI",
                   T0=p1.1$coefficients$estimate)
    eo <- c(paste0(" Fuzzy RD parameter  -0.2021",
                   "     0.1534      0.09753   (-0.5545, 0.1503)"),
            "Maximal leverage for fuzzy RD parameter: 0.003629",
            "First stage estimate: 0.3416 ")
    expect_equal(capture.output(print(p2, digits=4))[c(9, 14:15)],
                 eo)
})
