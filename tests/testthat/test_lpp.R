context("Test Inference at a point")

test_that("Inference at point agrees with RD", {
    d <- RDData(lee08, cutoff=0)
    rde <- RDHonest.fit(d, hp=5, M=2)
    dp <- LPPData(lee08[lee08$margin>0, ], point=0)
    dm <- LPPData(lee08[lee08$margin<0, ], point=0)
    pp <- LPPHonest.fit(dp, h=5, M=2)
    mm <- LPPHonest.fit(dm, h=5, M=2)

    expect_equal(pp$estimate-mm$estimate, rde$estimate)
    expect_equal(pp$sd^2+mm$sd^2, rde$sd^2)
    expect_equal(pp$maxbias+mm$maxbias, rde$maxbias)
    expect_equal(mm$eff.obs+pp$eff.obs, rde$eff.obs)
})

test_that("ROT bandwidth check", {
    ## Interior
    d <- LPPData(lee08, point=0)
    b1 <- ROTBW.fit(d, kern="uniform", order=1)
    b2 <- ROTBW.fit(d, kern="uniform", order=1, boundary=FALSE)
    expect_equal(b1, b2)
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
