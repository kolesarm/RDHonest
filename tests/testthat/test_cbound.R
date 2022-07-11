context("Test estimation of bound on C")

test_that("Test estimation of bound on M in Lee data", {
    r <- RDHonest(voteshare~margin, data=lee08)
    r1 <- RDSmoothnessBound(r, s=100, separate=TRUE, multiple=TRUE, sclass="T")
    expect_equal(r1$estimate, c(0.001491679, 0.000221391))
    r2 <- RDSmoothnessBound(r, s=100, separate=FALSE, multiple=TRUE, sclass="H")
    expect_equal(r2$conf.low, 0.03850396)
    r3 <- RDSmoothnessBound(r, s=100, separate=FALSE, multiple=FALSE,
                            sclass="H")
    expect_equal(r3$conf.low, 0.64174631)
    ## Should fail, s too big
    expect_error(lRDSmoothnessBound(rr, s=10000))
    ## Old implementation
    ## r2 <- capture.output(print(r, digits=5))
    ## expect_equal(c(r$po$Delta, r$ne$Delta),
    ##              c(0.002546658, 0.020467968))
    ## expect_identical(r2[c(5, 13)], c("Delta: 0.0025467, sd=0.0015372",
    ##                                  "Delta: 0.020468, sd=0.0073167"))
})
