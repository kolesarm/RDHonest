context("Test estimation of bound on C")

test_that("Test estimation of bound on M in Lee data", {

    d <- RDPrelimVar(RDData(lee08, cutoff=0), se.initial="IKEHW")
    r <- RDSmoothnessBound(d, s=100, separate=TRUE, multiple=TRUE, sclass="T")
    expect_equal(c(r$po$Delta, r$ne$Delta),
                 c(0.002546658, 0.020467968))
})
