context("Test Inference at a point")

test_that("Inference at point agrees with RD", {
    d <- RDData(lee08, cutoff=0)
    rde <- NPRHonest.fit(d, h=5, M=2)
    dp <- LPPData(lee08[lee08$margin>0, ], point=0)
    dm <- LPPData(lee08[lee08$margin<0, ], point=0)
    pp <- NPRHonest.fit(dp, h=5, M=2)
    mm <- NPRHonest.fit(dm, h=5, M=2)
    expect_equal(pp$estimate-mm$estimate, rde$estimate)
    expect_equal(pp$sd^2+mm$sd^2, rde$sd^2)
    expect_equal(pp$maxbias+mm$maxbias, rde$maxbias)
    expect_equal(mm$eff.obs+pp$eff.obs, rde$eff.obs)

    p2 <- LPPHonest(voteshare~margin, data=lee08, subset=margin>=0, h=5, M=2)
    ## 1:9 since something weird happens on codecov.io if lintr is included
    expect_equal(capture.output(print(pp))[1:9],
                 capture.output(print(p2))[6:14])
    ## Local constant yields infinite bias
    expect_equal(NPRHonest.fit(d, h=5, M=2, order=0)$maxbias, Inf)
    expect_equal(NPRHonest.fit(dp, h=10, M=0.1, order=0)$maxbias, Inf)

    ## Compare RDHonest and LPPHonest when order=2
    r <- NPRHonest.fit(d, h=7, M=2, order=2)
    rm <- NPRHonest.fit(dp, h=7, M=2, order=2)
    rp <- NPRHonest.fit(dm, h=7, M=2, order=2)
    expect_equal(r$maxbias, rm$maxbias+rp$maxbias)
    expect_equal(sqrt(rp$sd^2+rm$sd^2), r$sd)
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
    b1 <- ROTBW.fit(d, kern="uniform", order=1)
    b2 <- ROTBW.fit(d, kern="uniform", order=1, boundary=FALSE)
    expect_equal(b1, b2)
    expect_equal(ROTBW.fit(d, kern="triangular", order=1),
                 ROTBW.fit(d, kern=function(u) pmax(1-abs(u), 0), order=1))

    ## f0 using Silverman:
    f0 <- 0.0089934638
    C <- 9/8                            # nu0/(4*mu_2^2)
    ll <- lm(d$Y~d$X+I(d$X^2)+I(d$X^3)+I(d$X^4))
    h <- (C*sigma(ll)^2/(length(d$X)*f0*ll$coefficients[3]^2))^(1/5)
    expect_equal(b2, unname(h))

    dp <- LPPData(lee08[lee08$margin>0, ], point=0)
    bp1 <- ROTBW.fit(dp, kern="uniform", order=1)
    bp2 <- ROTBW.fit(dp, kern="uniform", order=1, boundary=TRUE)
    expect_equal(bp1, bp2)

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
    expect_equal(r[13], "Bandwidth: 13.41")
    r2 <- NPROptBW.fit(dp, M=2*Mh, opt.criterion="MSE")$h
    expect_identical(r2,  r1$h)
    expect_equal(r2,  13.4109132)
    expect_equal(unname(NPRHonest.fit(dp, M=Mh, kern="uniform",
                                      opt.criterion="FLCI")$upper),
                 55.24963853)

    ## Make sure we're getting positive worst-case bias
    leep <- lee08[lee08$margin>0, ]
    M <- NPR_MROT.fit(LPPData(leep, point = 20))
    r <- LPPHonest(voteshare ~ margin, data=leep, point=20, kern="uniform", M=M,
                   opt.criterion="MSE", sclass="H")
    expect_equal(r$maxbias, 0.2482525)
})
