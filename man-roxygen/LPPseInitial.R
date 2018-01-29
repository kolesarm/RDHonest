#' @param se.initial Method for estimating initial variance for computing
#'     optimal bandwidth. Ignored if data already contains estimate of variance.
#'
#' \describe{
#'
#' \item{"ROTEHW"}{Based on residuals from a local linear regression using a
#'              triangular kernel and ROT bandwidth}
#'
#' \item{"ROTdemeaned"}{Based on sum of squared deviations of outcome from
#'               estimate of intercept in local linear regression with
#'               triangular kenrel and ROT bandwidth}
#' }
