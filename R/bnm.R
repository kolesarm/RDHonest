#' Critical values for CIs based on a biased Gaussian estimator.
#'
#' Computes the critical value \eqn{cv_{1-alpha}(B)} that is needed to make the
#' confidence interval \eqn{X\pm cv} have coverage \eqn{1-alpha} if \eqn{X} is
#' Normally distributed with variance one and maximum bias at most \eqn{B}.
#'
#' @param B Maximum bias, vector of non-negative numbers.
#' @param alpha Determines CI level, \code{1-alpha}. Needs to be between 0 and
#' 1. Can be a vector of values.
#' @return Data frame with the following columns:
#' \describe{
#' \item{bias}{Value of bias as specfied by \code{bs}}
#' \item{alpha}{Value of \eqn{\alpha} as specified by \code{alpha}}
#' \item{cv}{Critical value}
#' \item{TeXDescription}{LaTeX-friendly description of current row}
#' }
#' @examples
#' # 90% critical value:
#' CVb(B = 1, alpha = 0.1)
#' CVb(B = c(0, 0.5, 1), alpha = c(0.05, 0.1))
#' @importFrom stats pnorm qnorm
#' @export
CVb <- function(B, alpha=0.05) {
    cv <- function(B, alpha = 0.05) {
        stopifnot(B >= 0 & alpha > 0 & alpha < 1)
        stats::uniroot(function(c) pnorm(c - B) - pnorm(- c - B) - (1-alpha),
                       c(1e-10, B - 2*qnorm(alpha / 2)),
                       tol=tol)$root
    }

    d <- data.frame(expand.grid(bias=B, alpha=alpha), cv=NA)
    d$cv <- sapply(1:nrow(d), function(j) cv(d$bias[j], d$alpha[j]))
    d$TeXDescription <- paste("$\\alpha=", d$alpha, "$", sep="")

    d

}
