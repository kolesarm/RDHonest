context("Test estimation of bound on C")

test_that("Test estimation of bound on M in Lee data", {
    r <- RDHonest(voteshare~margin, data=lee08)
    r1 <- RDSmoothnessBound(r, s=100, separate=TRUE, multiple=TRUE, sclass="T")
    expect_equal(r1$estimate, c(0.00147446, 0.00016410))
    r2 <- RDSmoothnessBound(r, s=100, separate=FALSE, multiple=TRUE, sclass="H")
    expect_equal(r2$conf.low, 0.017285255)
    r3 <- RDSmoothnessBound(r, s=100, separate=FALSE, multiple=FALSE,
                            sclass="H")
    expect_equal(r3$conf.low, 0.59509966)
    ## Should fail, s too big
    expect_error(RDSmoothnessBound(r, s=10000))
    r <- RDHonest(log(earnings)~yearat14, cutoff=1947, data=cghs, M=0.023,
                  h=3.4)
    expect_error(RDSmoothnessBound(r, s=20, separate=TRUE, multiple=FALSE))
    r4 <- RDSmoothnessBound(r, s=2, separate=TRUE, multiple=FALSE)
    expect_equal(r4$estimate, c(0.015359658, 0L))


    ## Old implementation
    ## r2 <- capture.output(print(r, digits=5))
    ## expect_equal(c(r$po$Delta, r$ne$Delta),
    ##              c(0.002546658, 0.020467968))
    ## expect_identical(r2[c(5, 13)], c("Delta: 0.0025467, sd=0.0015372",
    ##                                  "Delta: 0.020468, sd=0.0073167"))
})
