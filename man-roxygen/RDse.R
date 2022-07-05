#' @param se.method Vector with methods for estimating standard error of
#' estimate. If \code{NULL}, standard errors are not computed. The elements of
#' the vector can consist of the following methods:
#'
#' \describe{
#'     \item{"nn"}{Nearest neighbor method}
#'
#'     \item{"EHW"}{Eicker-Huber-White, with residuals from local regression
#'     (local polynomial estimators only).}
#'
#'    \item{"supplied.var"}{Use conditional variance supplied by \code{sigma2} or
#'         \code{d} instead of computing residuals}
#'
#' }
#' @param J Number of nearest neighbors, if "nn" is specified in
#'     \code{se.method}.
