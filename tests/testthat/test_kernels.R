context("Equivalent kernels")

test_that("Numeric and analytical equivalent kernels agree", {
    us <- matrix(c(-2, seq(-0.9, 0.9, by=0.05), 0, 2), nrow = 4)

    ## first  boundary
    Kt <- function(u) 2 * (1 - abs(u)) * (u <= 1) * (u >= 0)
    Ku <- function(u) (u <= 1) * (u >= 0)
    Ke <- function(u) (3 / 2) * (1 - u^2) * (u <= 1) * (u >= 0)

    for (j in 0:4) {
        expect_equal(EqKern(Ku, boundary = TRUE, order = j)(us),
                     EqKern("uniform", boundary = TRUE, order = j)(us))
        expect_equal(EqKern(Kt, boundary = TRUE, order = j)(us),
                     EqKern("triangular", boundary = TRUE, order = j)(us))
        expect_equal(EqKern(Ke, boundary = TRUE, order = j)(us),
                     EqKern("epanechnikov", boundary = TRUE, order = j)(us))
    }

    ## now interior
    Ku <- function(u) (u <= 1) * (u >= -1) / 2
    Kt <- function(u) (1 - abs(u)) * (u <= 1) * (u >= -1)
    Ke <- function(u) (3 / 4) * (1 - u^2) * (u <= 1) * (u >= -1)

    for (j in (0:4)) {
        expect_equal(EqKern(Ku, boundary = FALSE, order = j)(us),
                     EqKern("uniform", boundary = FALSE, order = j)(us))
        expect_equal(EqKern(Kt, boundary = FALSE, order = j)(us),
                     EqKern("triangular", boundary = FALSE, order = j)(us))
        expect_equal(EqKern(Ke, boundary = FALSE, order = j)(us),
                     EqKern("epanechnikov", boundary = FALSE, order = j)(us))
    }

})

test_that("Analytical kernel moments match numerical ones", {

    kernN <- cbind(kernC[, 1:3], kernC[, 4:13]*0)
    NumericalMoment <- function(moment, type) {
        NumericalMomentj <- function(j) {
            K <- EqKern(kernN$kernel[j], kernN$boundary[j],
                        kernN$order[j])
            KernMoment(K, moment, kernN$boundary[j], type = type)
        }

        vapply(seq_len(nrow(kernN)), NumericalMomentj, numeric(1))
    }

    for (j in 0:4) {
        kernN[[paste("mu", j, sep = "")]] <- NumericalMoment(j, "raw")
        kernN[[paste("nu", j, sep = "")]] <- NumericalMoment(j, "raw2")
        kernN[[paste("pi", j, sep = "")]] <- NumericalMoment(j, "absolute")
    }

    expect_equal(kernN[, 1:8], kernC[, 1:8]) # mu
    expect_equal(kernN[, 9:13], kernC[, 9:13]) # nu
    ## for pi, we allow some numerical imprecision
    expect_true(max(abs(kernN[, 14:18]- kernC[, 14:18]))<10e-6)
})
