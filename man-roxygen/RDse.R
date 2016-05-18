#' @param se.method Vector with methods for estimating standard error of
#'     estimate. If \code{NULL}, use estimate of sigma^2(x) supplied by
#'     \code{sigma2}. The elements of the vector can consist of the following
#'     method:
#'
#' \describe{
#'     \item{"nn"}{Nearest neighbor method}
#'     \item{"EHW"}{Eicker-Huber-White, with residuals from local regression.}
#'
#'     \item{"demeaned"}{Use EHW, but instead of using residuals, estimate
#'     \eqn{sigma^2_i} by subtracting the estimated intercept from the
#'     outcome (and not subtracting the estimated slope)}}
#' @param J Number of nearest neighbors, if "nn" is specified in \code{se.method}.
