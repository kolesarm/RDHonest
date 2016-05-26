context("Test bounded normal mean")

test_that("Sanity check for CVb", {

    alpha <- 0.07
    expect_equal(CVb(100, alpha=alpha)$cv, qnorm(1-alpha) + 100)
    expect_equal(CVb(0, alpha=alpha)$cv, qnorm(1-alpha/2))
})
