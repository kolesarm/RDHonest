#' Critical values for CIs based on a biased Gaussian estimator.
#'
#' Computes the critical value \eqn{cv_{1-\alpha}(B)} such that the confidence
#' interval \eqn{X\pm cv_{1-\alpha}(B)} will have coverage \eqn{1-\alpha}, where
#' \eqn{X} is Normally distributed with variance equal to \eqn{1} and maximum
#' bias at most \eqn{B}.
#'
#' @param B Maximum bias, vector of non-negative numbers.
#' @param alpha Determines CI level, \code{1-alpha}. Needs to be between 0 and
#'     1. Can be a vector of values.
#' @return Data frame with the following columns:
#'
#' \describe{
#'
#' \item{bias}{Value of bias as specified by \code{B}}
#'
#' \item{alpha}{Value of \eqn{\alpha} as specified by \code{alpha}}
#'
#' \item{cv}{Critical value}
#'
#' \item{TeXDescription}{LaTeX-friendly description of current row}}
#' @examples
#' ## 90% critical value:
#' CVb(B = 1, alpha = 0.1)
#' ## Returns data frame with 4 rows
#' CVb(B = c(0, 0.5, 1), alpha = c(0.05, 0.1))
#' @importFrom stats pnorm qnorm
#' @export
CVb <- function(B, alpha=0.05) {
    ## Single B
    cv <- function(B, alpha = 0.05) {
        if(is.na(B)) return(NA)
        stopifnot(B >= 0 & alpha > 0 & alpha < 1)
        sqrt(stats::qchisq(1-alpha, df = 1, ncp = B^2))
    }

    d <- data.frame(expand.grid(bias=B, alpha=alpha), cv=NA)
    d$cv <- vapply(seq_len(nrow(d)),
                   function(j) cv(d$bias[j], d$alpha[j]), numeric(1))
    d$TeXDescription <- paste("$\\alpha=", d$alpha, "$", sep="")

    d
}
