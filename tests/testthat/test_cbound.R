context("Test estimation of bound on C")

test_that("Test estimation of bound on M in Lee data", {

    d <- NPRPrelimVar.fit(RDData(lee08, cutoff=0), se.initial="EHW")
    r <- RDSmoothnessBound(d, s=100, separate=TRUE, multiple=TRUE, sclass="T")
    r2 <- capture.output(print(r, digits=5))
    expect_equal(c(r$po$Delta, r$ne$Delta),
                 c(0.002546658, 0.020467968))
    expect_identical(r2[c(5, 13)], c("Delta: 0.0025467, sd=0.0015372",
                                     "Delta: 0.020468, sd=0.0073167"))
})
