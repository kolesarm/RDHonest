context("Test Inference at a point")

test_that("Inference at point agrees with RD", {
    d <- NPRData(lee08, cutoff=0, "SRD")
    rde <- NPRHonest.fit(d, h=5, M=2)$coefficients
    dp <- NPRData(lee08[lee08$margin>=0, ], cutoff=0, "IP")
    dm <- NPRData(lee08[lee08$margin<0, ], cutoff=0, "IP")
    p0 <- NPRHonest.fit(dp, h=5, M=2)
    pp <- p0$coefficients
    mm <- NPRHonest.fit(dm, h=5, M=2)$coefficients
    expect_equal(pp$estimate-mm$estimate, rde$estimate)
    expect_equal(pp$std.error^2+mm$std.error^2, rde$std.error^2)
    expect_equal(pp$maximum.bias+mm$maximum.bias, rde$maximum.bias)



    p2 <- RDHonest(voteshare~margin, data=lee08, subset=margin>=0, h=5, M=2,
                   point.inference=TRUE)
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
    Mh <- RDHonest(mortHS~povrate, data=headst, cutoff=0, h=0)$coefficients$M
    expect_equal(Mh, 0.29939992)

    Mp <- RDHonest(mortHS~povrate, data=headst, subset=(povrate>=0),
                   cutoff=0, h=0, point.inference=TRUE)
    Mm <- RDHonest(mortHS~povrate, data=headst, subset=(povrate<0),
                   h=0, point.inference=TRUE)
    expect_equal(Mp$coefficients$M, Mh)
    expect_lt(Mm$coefficients$M, Mh)
})

test_that("ROT bandwidth check", {
    ## Interior
    d <- NPRData(lee08, cutoff=0, "IP")
    b1 <- ROTBW.fit(d, kern="uniform")

    ## f0 using Silverman:
    f0 <- 0.0089934638
    C <- 9/8                            # nu0/(4*mu_2^2)
    ll <- lm(d$Y~d$X+I(d$X^2)+I(d$X^3)+I(d$X^4))
    h <- (C*sigma(ll)^2/(length(d$X)*f0*ll$coefficients[3]^2))^(1/5)
    expect_equal(b1, unname(h))

    dp <- NPRData(lee08[lee08$margin>0, ], cutoff=0, "IP")
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
    rr <- RDHonest(voteshare ~ margin, data=lee08, subset=(margin>0),
                   kern="uniform", opt.criterion="FLCI", point.inference=TRUE)
    expect_equal(rr$coefficients$conf.high.onesided,
                 55.24963853)
    expect_equal(rr$coefficients$eff.obs, 858)
    expect_equal(rr$coefficients$eff.obs,
                 sum(lee08$margin<=rr$coefficients$bandwidth & lee08$margin>0))
    dp <- NPRData(lee08[lee08$margin>0, ], cutoff=0, "IP")
    Mh <- rr$coefficients$M

    re <- RDHonest(voteshare ~ margin, data=lee08, subset=(margin>0),
                   kern="uniform", opt.criterion="FLCI", point.inference=TRUE,
                   se.method="EHW")
    ## Should match regression
    rl <- lm(voteshare ~ margin, data=lee08,
             subset=margin>0 & margin<=rr$coefficients$bandwidth)
    XX <- model.matrix(rl)
    meat <- crossprod(XX, rl$residuals^2*XX)
    vl <- (solve(crossprod(XX)) %*% meat %*% solve(crossprod(XX)))[1, 1]
    expect_equal(sqrt(vl), re$coefficients$std.error)
    expect_equal(unname(rl$coefficients[1]), rr$coefficients$estimate)

    r1 <- RDHonest(voteshare~margin, data=lee08, subset=margin>0, M=2*Mh,
                   opt.criterion="MSE", point.inference=TRUE)
    r <- capture.output(print(r1, digits=4))
    expect_equal(r[11], "Bandwidth: 13.41, Kernel: triangular")
    r2 <- NPROptBW.fit(dp, M=2*Mh, opt.criterion="MSE")
    expect_identical(r2,  r1$coefficients$bandwidth)
    expect_lt(abs(r2- 13.4109133), 1e-6)

    ## Make sure we're getting positive worst-case bias
    r <- RDHonest(voteshare ~ margin, data=lee08, subset=(margin>0), cutoff=20,
                  kern="uniform", opt.criterion="MSE", point.inference=TRUE)
    expect_equal(r$coefficients$maximum.bias, 0.2482525)
})
