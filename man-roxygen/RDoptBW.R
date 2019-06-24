#' @param opt.criterion Optimality criterion that bandwidth is designed to
#'     optimize. The options are:
#'
#'    \describe{
#'
#'    \item{\code{"MSE"}}{Finite-sample maximum MSE}
#'
#'    \item{\code{"FLCI"}}{Length of (fixed-length) two-sided
#'        confidence intervals.}
#'
#'    \item{\code{"OCI"}}{Given quantile of excess length of one-sided
#'        confidence intervals}
#'
#'     }
#'
#'     The methods use conditional variance given by \code{sigma2}, if supplied.
#'     Otherwise, for the purpose of estimating the optimal bandwidth,
#'     conditional variance is estimated using the method specified by
#'     \code{se.initial}.
#' @param beta Determines quantile of excess length to optimize, if bandwidth
#'     optimizes given quantile of excess length of one-sided confidence
#'     intervals; otherwise ignored.
#' @param alpha determines confidence level, \code{1-alpha} for
#'     constructing/optimizing confidence intervals.
