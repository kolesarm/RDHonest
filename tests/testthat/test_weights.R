context("Test weighted RD")

test_that("Test weighting using cghs", {
    d <- cghs
    ## Make 10 groups
    d$mod <- floor(10*(d$earnings - floor(d$earnings)))
    ## Make cells by group and year
    d$cell <- d$mod/10+d$yearat14
    dd <- data.frame()
    for (j in unique(d$cell)) {
        dd <- rbind(dd, data.frame(y=mean(log(d$earnings)[d$cell==j]),
                                   x=mean(d$yearat14[d$cell==j]),
                                   weights=length(d$yearat14[d$cell==j])))
    }

    d1 <- RDData(dd, cutoff=1947)
    d2 <- RDData(data.frame(y=log(cghs$earnings), x=cghs$yearat14), cutoff=1947)
    ## Initial estimates
    r2 <- NPRreg.fit(d2, h=5, kern="triangular")
    r1 <- NPRreg.fit(d1, h=5, kern="triangular")

    ## Checks weights match
    wp1 <- vapply(unique(d1$Xp), function(j) sum(r1$wp[d1$Xp==j]), numeric(1))
    wp2 <- vapply(unique(d2$Xp), function(j) sum(r2$wp[d2$Xp==j]), numeric(1))
    wm1 <- vapply(unique(d1$Xm), function(j) sum(r1$wm[d1$Xm==j]), numeric(1))
    wm2 <- vapply(unique(d2$Xm), function(j) sum(r2$wm[d2$Xm==j]), numeric(1))
    expect_equal(wm1, wm2)
    expect_equal(wp1, wp2)
    ## Variance by hand
    np <- vapply(unique(d2$Xp), function(j) sum(d2$Xp==j), numeric(1))
    nm <- vapply(unique(d2$Xm), function(j) sum(d2$Xm==j), numeric(1))

    d2$sigma2p <- rep(mean(r2$sigma2p), length(d2$Xp))
    d2$sigma2m <- rep(mean(r2$sigma2m), length(d2$Xm))
    d1$sigma2p <- mean(r2$sigma2p)/d1$wp
    d1$sigma2m <- mean(r2$sigma2m)/d1$wm
    v1 <- sqrt(sum(wp1^2/np)*mean(r2$sigma2p)+sum(wm1^2/nm)*mean(r2$sigma2m))

    m2 <- NPRHonest.fit(d2, M=1, kern="triangular", h=5,
                        se.method="supplied.var")$coefficients
    m1 <- NPRHonest.fit(d1, M=1, kern="triangular", h=5,
                        se.method="supplied.var")$coefficients
    expect_equal(v1, m1$std.error)
    expect_equal(m1[2:6], m2[2:6])

    ## Same thing with RDHonest
    s1 <- RDHonest(log(earnings)~yearat14, cutoff=1947, h=5,
                   data=cghs, M=1)$coefficients
    s2 <- RDHonest(y~x, cutoff=1947, weights=weights, h=5,
                   data=dd, M=1)$coefficients
    expect_equal(c(s2$estimate, s2$maximum.bias),
                 c(s1$estimate, s1$maximum.bias))
    ## LPP Honest
    t1 <- RDHonest(log(earnings)~yearat14, cutoff=1947, h=5,
                   data=cghs[cghs$yearat14>=1947, ], M=1,
                   point.inference=TRUE)$coefficients
    t2 <- RDHonest(y~x, cutoff=1947, h=5, data=dd[dd$x>=1947, ], M=1,
                   weights=weights, point.inference=TRUE)$coefficients
    expect_equal(c(t2$estimate, t2$maximum.bias),
                 c(t1$estimate, t1$maximum.bias))
})
