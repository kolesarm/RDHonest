test_that("Test NN variance estimator", {
    r0 <- RDHonest(cn~retired | elig_year, data=rcp[1:5000, ],
                   M=c(1, 1), h=10)
    d <- r0$d

    s0 <- sigmaNN(d$X[d$p], d$Y[d$p, ], J=5)
    s1 <- sigmaNN(d$X[d$p], d$Y[d$p, 1], J=5)
    s2 <- sigmaNN(d$X[d$p], d$Y[d$p, 2], J=5)
    expect_equal(s0[, 1], s1)
    expect_equal(s0[, 4], s2)
})

test_that("Test LPreg", {
    r0 <- RDHonest(cn~retired | elig_year, data=rcp[1:5000, ],
                   M=c(1, 1), h=10)
    d <- r0$d
    d$sigma2 <- NA
    K <- EqKern("triangular", boundary=FALSE, order=0)
    r0e <- LPReg(d$X[d$m], d$Y[d$m, ], h=10, K, order=2, se.method="EHW")
    r0n <- LPReg(d$X[d$m], d$Y[d$m, ], h=10, K, order=2, se.method="nn", J=4)
    r1e <- LPReg(d$X[d$m], d$Y[d$m, 1], h=10, K, order=2, "EHW")
    r1n <- LPReg(d$X[d$m], d$Y[d$m, 1], h=10, K, order=2, "nn", J=4)
    r2e <- LPReg(d$X[d$m], d$Y[d$m, 2], h=10, K, order=2, se.method="EHW")
    r2n <- LPReg(d$X[d$m], d$Y[d$m, 2], h=10, K, order=2, se.method="nn", J=4)
    expect_lt(max(abs(r0e$theta- c(r1e$theta, r2e$theta))), 1e-10)
    expect_lt(max(abs(r0n$theta- c(r1n$theta, r2n$theta))), 1e-10)

    expect_equal(unname(r0e$sigma2[, c(1, 4)]),
                 cbind(r1e$sigma2, r2e$sigma2))
    expect_equal(unname(r0n$sigma2[, c(1, 4)]),
                 cbind(r1n$sigma2, r2n$sigma2))

    expect_equal(unname(r0e$var[c(1, 4)]), c(r1e$var, r2e$var))
    expect_equal(unname(r0n$var[c(1, 4)]), c(r1n$var, r2n$var))

    expect_equal(r0e$w, r1e$w)
    expect_equal(r0e$w, r1n$w)

    expect_identical(r0e$eff.obs, r1e$eff.obs)
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
    d <- RDHonest(log(cn)~retired|elig_year, data=rcp, M=c(1, 1), h=6)$d

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
