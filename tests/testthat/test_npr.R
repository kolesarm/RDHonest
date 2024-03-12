test_that("Test NN variance estimator", {
    r0 <- RDHonest(cn|retired~ elig_year, data=rcp[1:500, ],
                   M=c(1, 1), h=20)
    d <- r0$d

    s0 <- sigmaNN(d$X[d$p], d$Y[d$p, ], J=5)
    s1 <- sigmaNN(d$X[d$p], d$Y[d$p, 1], J=5)
    s2 <- sigmaNN(d$X[d$p], d$Y[d$p, 2], J=5)
    expect_equal(s0[, 1], s1)
    expect_equal(s0[, 4], s2)
})

test_that("Test LPreg", {
    d <- rcp[1:1000, ]
    r0 <- RDHonest(cn ~ elig_year, data=d, M=1, h=10, J=4)
    r0m <- RDHonest(cn ~ elig_year, data=d, M=1, h=10, subset=elig_year<0,
                    point.inference=TRUE, J=4)
    r0p <- RDHonest(cn ~ elig_year, data=d, M=1, h=10, subset=elig_year>=0,
                    point.inference=TRUE, J=4)
    expect_equal(r0p$coefficients$estimate-r0m$coefficients$estimate,
                 r0$coefficients$estimate)
    expect_equal(range(c(r0m$data$sigma2,r0p$data$sigma2)-r0$data$sigma2),
                 c(0, 0))
    expect_equal(as.numeric(r0m$coefficients[c(2:10)]),
                 c(20139.543667707,  2097.006669471,     6.641178638,
                   16029.465508831, 24249.621826583, 16683.633463048,
                   23595.453872367,    10L,   125.503763304))
    NPReg(r0$data, h, kern="triangular", order=2, se.method="nn", J=3)
    rr <- NPReg(r0$data, h, kern="triangular", order=2, se.method="nn", J=6)
    expect_equal(as.numeric(rr[1:2]),  c(-5395.040191, 9754.921083))
})

test_that("Test NPReg", {
    ## Replicate Ludwig-Miller
    lumi <- headst[!is.na(headst$mortHS), c("mortHS", "povrate")]
    mort <- RDHonest(mortHS~povrate, data=lumi, M=1)$d
    mortm <- RDHonest(mortHS~povrate, data=lumi, M=1,
                      subset=povrate<0, point.inference=TRUE)$d
    mortp <- RDHonest(mortHS~povrate, data=lumi, M=1,
                      subset=povrate>=0, point.inference=TRUE)$d
    t1 <- data.frame()
    t2 <- data.frame()
    bws <- c(9, 18, 36)
    for (bw in bws) {
        rfe <- NPReg(mort, bw, "uniform", order=1, se.method="EHW")
        rfn <- NPReg(mort, bw, "uniform", order=1, se.method="nn")
        t1 <- rbind(t1, data.frame(bw=bw, estimate=rfe$estimate,
                                   EHW=rfe$se, nn=rfn$se))
        rme <- NPReg(mortm, bw, "uniform", order=1, se.method="EHW")
        rmn <- NPReg(mortm, bw, "uniform", order=1, se.method="nn")
        rpe <- NPReg(mortp, bw, "uniform", order=1, se.method="EHW")
        rpn <- NPReg(mortp, bw, "uniform", order=1, se.method="nn")
        t2 <- rbind(t2, data.frame(bw=bw, estimate=rpe$estimate-rme$estimate,
                                   EHW=sqrt(rpe$se^2+rme$se^2),
                                   nn=sqrt(rpn$se^2+rmn$se^2)))
    }
    exp <- rbind(c(9, -1.895234221, 0.9801414664, 1.0381954004),
                 c(18, -1.198258131, 0.6608206598, 0.6955768728),
                 c(36, -1.113938949, 0.5001394599, 0.5223659480))
    expect_equal(unname(as.matrix(t1)), exp)
    expect_lt(max(abs(t1-t2)), 1e-11)

    ## Replicate Battistin et al
    df <- RDHonest(retired~elig_year, data=rcp, M=1, h=6)$d
    dr <- RDHonest(log(cn)~elig_year, data=rcp, M=1, h=6)$d
    d <- RDHonest(log(cn)|retired~elig_year, data=rcp, M=c(1, 1), h=6)$d

    rf <- NPReg(df, 10, "uniform", order=1, se.method="EHW")
    rr <- NPReg(dr, 10, "uniform", order=1, se.method="EHW")
    expect_equal(unname(c(rf$estimate, round(rf$se, 4))),
                 c(0.43148435, 0.0181))
    expect_equal(unname(c(rr$estimate, round(rr$se, 4))),
                 c(-0.03550599496, 0.0211))
    r1 <- NPReg(d, 10, "uniform", order=1, se.method="EHW")
    expect_equal(unname(c(r1$estimate, r1$se)),
                 c(-0.08228802393, 0.04830389419))
    ## Numbers from
    ## dt <- cbind(logcn=log(rcp[, 6]), rcp[, c(3, 2)], Z=rcp$elig_year>=0)
    ## r4 <- AER::ivreg(logcn ~ retired + Z:elig_year + elig_year |
    ##                      Z:elig_year + elig_year + Z,
    ##                  subset=(abs(elig_year)<=10), data=dt)
    ## summary(r4, vcov = sandwich::sandwich,
    ##         diagnostics=TRUE)$coefficients[2, 1:2]
    expect_lt(abs(r1$estimate- rr$estimate/rf$estimate), 1e-10)
    r2 <- NPReg(d, 5, function(x) abs(x)<=1, order=0, se.method="EHW")
    ## r3 <- AER::ivreg(logcn~retired | Z, subset=(abs(elig_year)<=5), data=dt)
    ## summary(r3, vcov = sandwich::sandwich,
    ## diagnostics = TRUE)$coefficients[2, 1:2]
    expect_equal(unname(c(r2$estimate, r2$se)),
                 c(-0.14157767658, 0.02562096397))
})
