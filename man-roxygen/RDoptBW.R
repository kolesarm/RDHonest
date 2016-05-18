#' @param opt.criterion Optimality criterion that bandwidth is designed to
#'     optimize. It can either be based on exact finite-sample maximum bias and
#'     finite-sample estimate of variance, or asymptotic approximations to the
#'     bias and variance. The options are:
#'
#'    \describe{
#'
#'    \item{\code{"MSE"}, \code{"aMSE"}}{Finite-sample or asymtotic maximum
#'        MSE}
#'
#'    \item{\code{"FLCI"}, \code{"aFLCI"}}{Length of (fixed-length) two-sided
#'        confidence intervals.}
#'
#'    \item{\code{"OCI"}, \code{"aOCI"}}{Given quantile of excess length of
#'        one-sided confidence intervals}
#'
#'     }
#'
#'     The finite-sample methods use conditional variance given by
#'     \code{sigma2}, if supplied. Otherwise, for the purpose of estimating the
#'     optimal bandwidth, conditional variance is assumed homoscedastic, and
#'     estimated using a nearest neighbor estimator.
#' @param bw.equal logical specifying whether bandwidths on either side of
#'     cutoff should be constrainted to equal to each other.
#' @param beta Determines quantile of excess length to optimize, if bandwidth
#'     optimizes given quantile of excess length of one-sided confidence
#'     intervals.
#' @param alpha determines confidence level, \code{1-alpha} for
#'     constructing/optimizing confidence intervals.
