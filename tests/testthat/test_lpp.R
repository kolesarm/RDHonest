test_that("Inference at point agrees with RD", {
    lees <- lee08[(1:1000)*6, ]
    rde <- RDHonest(voteshare~margin, data=lees, M=2, h=10)
    pp <- RDHonest(voteshare~margin, data=lees, M=2, h=10,
                   subset=I(margin>=0), point.inference=TRUE)
    mm <- RDHonest(voteshare~margin, data=lees, M=2, h=10,
                   subset=I(margin<0), point.inference=TRUE)
    expect_equal(pp$coefficients[2]-mm$coefficients[2], rde$coefficients[2])
    expect_equal(pp$coefficients[3]^2+mm$coefficients[3]^2,
                 rde$coefficients[3]^2)
    expect_equal(pp$coefficients[4]+mm$coefficients[4], rde$coefficients[4])

    r <- NPRHonest(rde$d, h=7, M=2)$coefficients
    rm <- NPRHonest(mm$d, h=7, M=2)$coefficients
    rp <- NPRHonest(pp$d, h=7, M=2)$coefficients
    expect_equal(r$maximum.bias,
                 rm$maximum.bias+rp$maximum.bias)
    expect_equal(sqrt(rp$std.error^2+rm$std.error^2), r$std.error)
})

test_that("MROT matches paper", {
    suppressMessages(Mh <- RDHonest(mortHS~povrate, data=headst,
                                    cutoff=0, h=0)$coefficients$M)
    expect_equal(Mh, 0.29939992)

    suppressMessages(Mp <- RDHonest(mortHS~povrate, data=headst,
                                    subset = (povrate>=0), cutoff=0, h=0,
                                    point.inference=TRUE))
    suppressMessages(Mm <- RDHonest(mortHS~povrate, data=headst,
                                    subset = (povrate<0), h=0,
                                    point.inference=TRUE))
    expect_equal(Mp$coefficients$M, Mh)
    expect_lt(Mm$coefficients$M, Mh)
})

test_that("ROT bandwidth check", {
    ## Interior
    r0 <- RDHonest(voteshare~margin, data=lee08, point.inference=TRUE, M=10)
    d <- r0$d
    b1 <- ROTBW(d, kern="uniform")

    ## f0 using Silverman:
    f0 <- 0.0089934638
    C <- 9/8                            # nu0/(4*mu_2^2)
    ll <- lm(d$Y~d$X+I(d$X^2)+I(d$X^3)+I(d$X^4))
    h <- (C*sigma(ll)^2 / (length(d$X)*f0*ll$coefficients[3]^2))^(1/5)
    expect_equal(b1, unname(h))

    r1 <- RDHonest(voteshare~margin, data=lee08, point.inference=TRUE,
                   subset=margin>0, M=1)$d
    bp1 <- ROTBW(r1, kern="uniform")

    ## f0 using Silverman:
    f0 <- 0.0079735105
    C <- 36                            # nu0/(4*mu_2^2)
    ll <- lm(r1$Y~r1$X+I(r1$X^2)+I(r1$X^3)+I(r1$X^4))
    der <- unname(ll$coefficients[3])
    h <- (C*sigma(ll)^2 / (length(r1$X)*f0*der^2))^(1/5)
    expect_equal(bp1, h)
})

test_that("Optimal bandwidth calculations", {
    expect_message(rr <- RDHonest(voteshare ~ margin, data=lee08,
                                  subset = (margin>0),
                                  kern="uniform", opt.criterion="FLCI",
                                  point.inference=TRUE))
    expect_equal(rr$coefficients$conf.high.onesided,
                 55.24963853)
    expect_equal(rr$coefficients$eff.obs, 858)
    expect_equal(rr$coefficients$eff.obs,
                 sum(lee08$margin<=rr$coefficients$bandwidth & lee08$margin>0))
    Mh <- rr$coefficients$M

    expect_message(re <- RDHonest(voteshare ~ margin, data=lee08,
                                  subset = (margin>0),
                                  kern="uniform", opt.criterion="FLCI",
                                  point.inference=TRUE, se.method="EHW"))
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

    re$d$sigma2 <- NULL
    r2 <- OptBW(re$d, M=2*Mh, opt.criterion="MSE")
    expect_identical(r2,  r1$coefficients$bandwidth)
    expect_lt(abs(r2- 13.4109133), 1e-6)

    ## Make sure we're getting positive worst-case bias
    expect_message(r <- RDHonest(voteshare ~ margin, data=lee08,
                                 subset = (margin>0), cutoff=20, kern="uniform",
                                 opt.criterion="MSE", point.inference=TRUE))
    expect_equal(r$coefficients$maximum.bias, 0.248252496)
})
