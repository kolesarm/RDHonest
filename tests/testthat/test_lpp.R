context("Test Inference at a point")

test_that("Inference at point agrees with RD", {
    d <- RDData(lee08, cutoff=0)
    rde <- NPRHonest.fit(d, h=5, M=2)$coefficients
    dp <- LPPData(lee08[lee08$margin>0, ], point=0)
    dm <- LPPData(lee08[lee08$margin<0, ], point=0)
    p0 <- NPRHonest.fit(dp, h=5, M=2)
    pp <- p0$coefficients
    mm <- NPRHonest.fit(dm, h=5, M=2)$coefficients
    expect_equal(pp$estimate-mm$estimate, rde$estimate)
    expect_equal(pp$std.error^2+mm$std.error^2, rde$std.error^2)
    expect_equal(pp$maximum.bias+mm$maximum.bias, rde$maximum.bias)
    expect_equal(mm$eff.obs+pp$eff.obs, rde$eff.obs)

    p2 <- LPPHonest(voteshare~margin, data=lee08, subset=margin>=0, h=5, M=2)
    ## 1:9 since something weird happens on codecov.io if lintr is included
    expect_equal(capture.output(print(p0))[1:7],
                 capture.output(print(p2))[6:12])

    r <- NPRHonest.fit(d, h=7, M=2)$coefficients
    rm <- NPRHonest.fit(dp, h=7, M=2)$coefficients
    rp <- NPRHonest.fit(dm, h=7, M=2)$coefficients
    expect_equal(r$maximum.bias,
                 rm$maximum.bias+rp$maximum.bias)
    expect_equal(sqrt(rp$std.error^2+rm$std.error^2), r$std.error)
    ## Silverman variance should approximately match
})

test_that("MROT matches paper", {
    ## HS Mortality
    lumi <- headst[!is.na(headst$mortHS), c("mortHS", "povrate60")]
    mort <- RDData(lumi, cutoff=0)
    Mh <- NPR_MROT.fit(mort)
    dp <- LPPData(lumi[lumi$povrate>=0, ], point=0)
    dm <- LPPData(lumi[lumi$povrate<0, ], point=0)

    expect_equal(Mh, 0.29939992)
    expect_equal(Mh, max(NPR_MROT.fit(dp), NPR_MROT.fit(dm)))
})

test_that("ROT bandwidth check", {
    ## Interior
    d <- LPPData(lee08, point=0)
    b1 <- ROTBW.fit(d, kern="uniform")
    expect_equal(ROTBW.fit(d, kern="triangular"),
                 ROTBW.fit(d, kern=function(u) pmax(1-abs(u), 0)))

    ## f0 using Silverman:
    f0 <- 0.0089934638
    C <- 9/8                            # nu0/(4*mu_2^2)
    ll <- lm(d$Y~d$X+I(d$X^2)+I(d$X^3)+I(d$X^4))
    h <- (C*sigma(ll)^2/(length(d$X)*f0*ll$coefficients[3]^2))^(1/5)
    expect_equal(b1, unname(h))

    dp <- LPPData(lee08[lee08$margin>0, ], point=0)
    bp1 <- ROTBW.fit(dp, kern="uniform")

    ## f0 using Silverman:
    f0 <- 0.0079735105
    C <- 36                            # nu0/(4*mu_2^2)
    ll <- lm(dp$Y~dp$X+I(dp$X^2)+I(dp$X^3)+I(dp$X^4))
    der <- unname(ll$coefficients[3])
    h <- (C*sigma(ll)^2/(length(dp$X)*f0*der^2))^(1/5)
    expect_equal(bp1, h)
})

test_that("Optimal bandwidth calculations", {
    dp <- LPPData(lee08[lee08$margin>0, ], point=0)
    Mh <- NPR_MROT.fit(dp)
    r1 <- LPPHonest(voteshare~margin, data=lee08, subset=margin>0, point=0,
                    M=2*Mh, opt.criterion="MSE")
    r <- capture.output(print(r1, digits=4))
    expect_equal(r[11], "Bandwidth: 13.41")
    r2 <- NPROptBW.fit(dp, M=2*Mh, opt.criterion="MSE")$h
    expect_identical(r2,  r1$coefficients$bandwidth)
    expect_lt(abs(r2- 13.4109133), 1e-6)
    rr <- NPRHonest.fit(dp, M=Mh, kern="uniform", opt.criterion="FLCI")
    expect_equal(rr$coefficients$conf.high.onesided,
                 55.24963853)

    ## Make sure we're getting positive worst-case bias
    leep <- lee08[lee08$margin>0, ]
    M <- NPR_MROT.fit(LPPData(leep, point = 20))
    r <- LPPHonest(voteshare ~ margin, data=leep, point=20, kern="uniform", M=M,
                   opt.criterion="MSE", sclass="H")
    expect_equal(r$coefficients$maximum.bias, 0.2482525)
})
