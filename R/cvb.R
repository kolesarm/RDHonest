#' Critical values for CIs based on a biased Gaussian estimator.
#'
#' Computes the critical value \eqn{cv_{1-\alpha}(B)}{cv_{1-alpha}(B)} such that
#' the confidence interval \eqn{X\pm cv_{1-\alpha}(B)}{X +/- cv_{1-alpha}(B)}
#' has coverage \eqn{1-\alpha}{1-alpha}, where the estimator \eqn{X} is normally
#' distributed with variance equal to \eqn{1} and maximum bias at most \eqn{B}.
#'
#' @param B Maximum bias, vector of non-negative numbers.
#' @param alpha Determines CI level, \eqn{1-\alpha}{1-alpha}. Scalar between 0
#'     and 1.
#' @return Vector of critical values, one for each value of maximum bias
#'     supplied by \code{B}.
#' @examples
#' ## 90% critical value:
#' CVb(B = 1, alpha = 0.1)
#' ## Usual 95% critical value
#' CVb(0)
#' ## Returns vector with 3 critical values
#' CVb(B = c(0, 0.5, 1), alpha = 0.05)
#' @export
CVb <- function(B, alpha=0.05) {
    cv <- B
    ## Take care of missing values
    cv <- function(B, alpha) {
        if (is.na(B)) return(NA)
        stopifnot(B >= 0 && alpha > 0 && alpha < 1)
        if (B<10)
            return(sqrt(stats::qchisq(1-alpha, df = 1, ncp = B^2)))
        else
            return(B+stats::qnorm(1-alpha))
    }
    vapply(seq_along(B), function(j) cv(B[j], alpha), numeric(1))
}
