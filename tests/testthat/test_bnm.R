context("Test bounded normal mean")

test_that("Sanity check for CVb for big and small values", {

    alpha <- 0.07
    expect_equal(CVb(100, alpha=alpha)$cv, qnorm(1-alpha) + 100)
    expect_equal(CVb(10000, alpha=alpha)$cv, qnorm(1-alpha) + 10000)
    expect_error(CVb(10000, alpha=0)$cv)
    expect_equal(CVb(0, alpha=alpha)$cv, qnorm(1-alpha/2))
})

test_that("CVb with NA values", {

    alpha <- c(0.07, 0.17)
    expect_equal(CVb(c(100, NA), alpha=alpha)$cv,
                     c(qnorm(1-alpha[1]) + 100, NA,
                       qnorm(1-alpha[2]) + 100, NA))
})
