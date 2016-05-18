#' @param se.method Vector with methods for estimating standard error of
#' estimate. If \code{NULL}, standard errors are not computed. The elements of
#' the vector can consist of the following methods:
#'
#' \describe{
#'     \item{"nn"}{Nearest neighbor method}
#'     \item{"EHW"}{Eicker-Huber-White, with residuals from local regression.}
#'     \item{"demeaned"}{Use EHW, but instead of using residuals, estimate
#'         \eqn{sigma^2_i} by subtracting the estimated intercept from the
#'         outcome (and not subtracting the estimated slope)}
#'     \item{"plugin"}{Plug-in estimate based on asymptotic variance.}
#'     \item{"supplied.var"}{Use EHW with conditional variance supplied by
#'         \code{sigma2} / \code{d} instead of computing residuals}
#' }
#' @param J Number of nearest neighbors, if "nn" is specified in \code{se.method}.
