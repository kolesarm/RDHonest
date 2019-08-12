context("Test NPR")

test_that("Test NN variance estimator", {
    d <- FRDData(rcp[1:5000, c(6, 3, 2)], cutoff=0)
    s0 <- sigmaNN(d$Xp, d$Yp, J=5)
    s1 <- sigmaNN(d$Xp, d$Yp[, 1], J=5)
    s2 <- sigmaNN(d$Xp, d$Yp[, 2], J=5)
    expect_equal(s0[, 1], s1)
    expect_equal(s0[, 4], s2)

})

test_that("Test LPreg", {
    d <- FRDData(rcp[1:5000, c(6, 3, 2)], cutoff=0)
    K <- EqKern("triangular", boundary=FALSE, order=0)
    r0 <- LPReg(d$Xm, d$Ym, h=10, K, order=2, se.method=c("EHW", "nn"), J=4)
    r1 <- LPReg(d$Xm, d$Ym[, 1], h=10, K, order=2, se.method=c("EHW", "nn"),
                J=4)
    r2 <- LPReg(d$Xm, d$Ym[, 2], h=10, K, order=2, se.method=c("EHW", "nn"),
                J=4)
    expect_identical(r0$theta, c(r1$theta, r2$theta))
    expect_equal(r0$sigma2[, c(1, 4)], cbind(r1$sigma2, r2$sigma2))
    expect_identical(r0$var[c(1, 4), ], rbind(r1$var, r2$var))
    expect_equal(r0$w, r1$w)
    expect_identical(r0$eff.obs, r1$eff.obs)
})

test_that("Test NPRreg", {
    ## Replicate Ludwig-Miller
    lumi <- headst[!is.na(headst$mortHS), c("mortHS", "povrate60")]
    mort <- RDData(lumi, cutoff=0)
    mortm <- LPPData(lumi[lumi$povrate<0, ], point=0)
    mortp <- LPPData(lumi[lumi$povrate>=0, ], point=0)
    t1 <- data.frame()
    t2 <- data.frame()
    bws <- c(9, 18, 36)
    for (bw in bws) {
        rf <- NPRreg.fit(mort, h=bw, kern="uniform",
                         order=1, se.method=c("EHW", "nn"), J=3)
        se <- unname(rf$se[c("EHW", "nn")])
        t1 <- rbind(t1, data.frame(bw=bw, estimate=rf$estimate,
                                   EHW=se[1], nn=se[2]))

        rm <- NPRreg.fit(mortm, h=bw, kern="uniform",
                         order=1, se.method=c("EHW", "nn"), J=3)
        rp <- NPRreg.fit(mortp, h=bw, kern="uniform",
                         order=1, se.method=c("EHW", "nn"), J=3)
        se <- unname(sqrt(rm$se[c("EHW", "nn")]^2 + rp$se[c("EHW", "nn")]^2))
        t2 <- rbind(t2, data.frame(bw=bw, estimate=rp$estimate-rm$estimate,
                                   EHW=se[1], nn=se[2]))

    }
    exp <- rbind(c(9, -1.895234221, 0.9801414664, 1.0381954004),
                 c(18, -1.198258131, 0.6608206598, 0.6955768728),
                 c(36, -1.113938949, 0.5001394599, 0.5223659480))
    expect_equal(unname(as.matrix(t1)), exp)
    expect_identical(max(abs(round(t1-t2, 14))), 0)
    ## Replicate Battistin et al
    df <- RDData(cbind(rcp[, c(3, 2)]), cutoff=0)
    dr <- RDData(cbind(logcn=log(rcp[, 6]), rcp[, 2, drop=FALSE]), cutoff=0)
    d <- FRDData(cbind(logcn=log(rcp[, 6]), rcp[, c(3, 2)]), cutoff=0)
    rf <- NPRreg.fit(df, h=10, kern="uniform", order=1, se.method="EHW")
    rr <- NPRreg.fit(dr, h=10, kern="uniform", order=1, se.method="EHW")
    expect_equal(unname(c(rf$estimate, round(rf$se["EHW"], 4))),
                 c(0.43148435, 0.0181))
    expect_equal(unname(c(rr$estimate, round(rr$se["EHW"], 4))),
                 c(-0.0355060, 0.0211))
    r1 <- NPRreg.fit(d, h=10, kern="uniform", order=1, se.method="EHW")
    expect_equal(c(r1$estimate, unname(r1$se["EHW"])),
                 c(-0.08228802, 0.0483039))
    ## Numbers from
    ## dt <- cbind(logcn=log(rcp[, 6]), rcp[, c(3, 2)], Z=rcp$elig_year>=0)
    ## r4 <- AER::ivreg(logcn ~ retired + Z:elig_year + elig_year |
    ##                      Z:elig_year + elig_year + Z,
    ##                  subset=(abs(elig_year)<=10), data=dt)
    ## summary(r4, vcov = sandwich::sandwich,
    ##         diagnostics=TRUE)$coefficients[2, 1:2]
    expect_identical(r1$estimate, rr$estimate/rf$estimate)
    r2 <- NPRreg.fit(d, h=5, kern=function(x) abs(x)<=1, order=0,
                     se.method="EHW")
    ## r3 <- AER::ivreg(logcn~retired | Z, subset=(abs(elig_year)<=5), data=dt)
    ## summary(r3, vcov = sandwich::sandwich,
    ## diagnostics = TRUE)$coefficients[2, 1:2]
    expect_equal(c(r2$estimate, unname(r2$se["EHW"])),
                 c(-0.14157767, 0.02562096))


    ## bw too narrow
    expect_warning(NPRreg.fit(d, h=2, kern="uniform", order=2))
    expect_warning(NPRreg.fit(RDData(lee08, cutoff=0), h=0.1, order=2))
})
