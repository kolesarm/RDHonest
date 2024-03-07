test_that("Test clustering formulas", {
    ## Everyone in their own cluster
    expect_message(s0 <- RDHonest(voteshare~margin, data=lee08,
                                  se.method="EHW"))
    expect_message(s0c <- RDHonest(voteshare~margin, data=lee08,
                                   se.method="EHW",
                                   clusterid=seq_along(lee08$margin)))
    expect_equal(s0$coefficients, s0c$coefficients)
    rr <- rcp[1:1000, ]
    expect_message(f0 <- RDHonest(c~retired | elig_year, data=rr,
                                  se.method="EHW"))
    expect_message(f0c <- RDHonest(c~retired | elig_year, data=rr,
                                   clusterid=seq_along(rr$c)))
    expect_equal(f0$coefficients, f0c$coefficients)
    expect_message(p0 <- RDHonest(c~elig_year, data=rr,
                                  se.method="EHW", point.inference=TRUE))
    expect_message(p0c <- RDHonest(c~elig_year, data=rr,
                                   clusterid=seq_along(rr$c),
                                   point.inference=TRUE))
    expect_equal(p0$coefficients, p0c$coefficients)

    ## Check we match sandwich::vcovCL: uniform + triangular
    set.seed(42)
    clusterid <- sample(1:50, NROW(lee08), replace=TRUE)
    expect_message(s1h <- RDHonest(clusterid~margin, data=lee08,
                                   se.method="EHW", h=5, kern="uniform"))
    expect_message(s1c <- RDHonest(clusterid~margin, data=lee08,
                                   se.method="EHW", h=5,
                                   kern="uniform", clusterid=clusterid))
    m1 <- lm(clusterid~margin*I(margin>0), data=lee08, subset=abs(margin)<=5)
    ## Should equal
    ## sqrt(sum(s1c$data$est_w[s1c$data$est_w!=0]^2*m1$residuals^2))
    ## or sqrt(sandwich::vcovHC(m1, type="HC0")[3, 3])
    expect_equal(as.numeric(s1h$coefficients[2:3]), c(-7.83271276, 2.388004382))
    ## sandwich::vcovCL(m1, type="HC0", cluster=clusterid[abs(lee08$margin)<=5],
    ##                  cadjust=FALSE)[3, 3]
    expect_equal(as.numeric(s1c$coefficients[3])^2, 6.523026427)
    expect_message(s2c <- RDHonest(clusterid~margin, data=lee08,
                                   se.method="EHW", clusterid=clusterid))
    h <- s2c$coefficients$bandwidth
    m2 <- lm(clusterid~margin*I(margin>0), data=lee08,
             subset=abs(margin)<=h,
             weights = (1 - abs(margin/h)))
    ## sqrt(sandwich::vcovCL(m2, type="HC0",
    ##                       cluster=clusterid[abs(lee08$margin)<=h],
    ##                       cadjust=FALSE)[3, 3])
    expect_equal(as.numeric(s2c$coefficients[2:3]),
                 c(unname(m2$coefficients[3]), 1.64034744815))

    set.seed(42)
    clusterid <- sample(1:20, NROW(rr), replace=TRUE)
    expect_message(p1c <- RDHonest(c~elig_year, data=rr, clusterid=clusterid,
                                   point.inference=TRUE, kern="uniform"))
    expect_message(p1h <- RDHonest(c~elig_year, data=rr, se.method="EHW",
                                   point.inference=TRUE, kern="uniform"))
    m2 <- lm(c~elig_year, data=rr, subset=abs(elig_year)<=5)
    expect_equal(p1c$coefficients$estimate, unname(m2$coefficients[1]))
    ## sqrt(sandwich::vcovHC(m2, type="HC0")[1, 1]) 792.2299439
    ## sandwich::vcovCL(m2, type="HC0", cluster=clusterid[abs(rr$elig_year)<=5],
    ##                       cadjust=FALSE)[1, 1] 651509.586
    expect_equal(p1c$coefficients$std.error^2, 723178.2655)
    expect_equal(p1h$coefficients$std.error, 792.2299439)

    expect_message(p2c <- RDHonest(c~elig_year, data=rr, clusterid=clusterid,
                                   point.inference=TRUE, h=8.27778063886))
    h <- p2c$coefficients$bandwidth
    m3 <- lm(c~elig_year, data=rr, subset=abs(elig_year)<=h,
             weights = (1 - abs(elig_year/h)))
    ## sqrt(sandwich::vcovCL(m3, type="HC0",
    ##                       cluster=clusterid[abs(rr$elig_year)<=h],
    ##                       cadjust=FALSE)[1, 1])
    expect_equal(as.numeric(p2c$coefficients[2:3]),
                 c(unname(m3$coefficients[1]), 919.162944))

    expect_message(f1c <- RDHonest(c~retired | elig_year, data=rr,
                                   clusterid=clusterid,
                                   kern="uniform"))
    expect_message(f2c <- RDHonest(cn~retired | elig_year, data=rr,
                                   clusterid=clusterid,
                                   kern="uniform", h=4))
    ## m3 <- AER::ivreg(c~retired+elig_year+elig_year:I(elig_year>0) |
    ##                      elig_year*I(elig_year>0),
    ##                  data=rr, subset=abs(elig_year)<=5)
    ## m3$coefficients[2] # -18532.2583128
    ## sqrt(sandwich::vcovCL(m3, type="HC0",
    ##                       cluster=clusterid[abs(rr$elig_year)<=5],
    ##                       cadjust=FALSE)[2,2]) # 13709.23278
    expect_equal(as.numeric(f1c$coefficients[2:3]),
                 c(-18532.2583128, 13709.23278))
    expect_equal(as.numeric(f2c$coefficients[2:3]),
                 c(-13623.08923, 14286.50171))
    ## Triangular kernel doesn't work with sandwich + ivreg

    ## Check preliminary bw calculations and rho
    r0 <- RDHonest(log(c)~elig_year, h=Inf, M=1, data=rr, clusterid=clusterid)
    expect_equal(PrelimVar(r0$d)$rho, -0.0008690793197)
})
