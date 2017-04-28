context("Test estimation of bound on C")

test_that("Test estimation of bound on M in Lee data", {

    d <- RDprelimVar(RDData(lee08, cutoff=0), se.initial="IKEHW")

    expect_equal(unname(RDMbound(d, 100)/2),
                 c(0.00635928,  0.00652546,  0.00568134,
                   0.00301211, -0.00342979,  0.00377666))
})
