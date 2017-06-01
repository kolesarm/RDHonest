#' @param se.initial Method for estimating initial variance for computing
#'     optimal bandwidth. Ignored data already contains estimate of variance.
#'
#' \describe{
#'
#' \item{"IKEHW"}{Based on residuals from a local linear regression using a
#'               triangular kenrel and IK bandwidth}
#'
#' \item{"IKdemeaned"}{Based on sum of squared deviations of outcome from
#'               estimate of intercept in local linear regression with
#'               triangular kenrel and IK bandwidth}
#'
#' \item{"Silverman"}{Use residuals from local constant regression with uniform
#' kernel and bandwidth selected using Silverman's rule of thumb, as in Equation
#' (14) in IK}

#' \item{"SilvermanNN"}{Use nearest neighbor estimates, rather than residuals}
#'
#' }