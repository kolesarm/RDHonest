#' @param se.initial Method for estimating initial variance for computing
#'     optimal bandwidth. Ignored if data already contains estimate of variance.
#'
#' \describe{
#'
#' \item{"IKEHW"}{Based on residuals from a local linear regression using a
#'               triangular kernel, with bandwidth corresponding to the Imbens
#'               and Kalyanaraman bandwidth for the reduced form}
#'
#' \item{"NN"}{Use nearest neighbor estimates, without assuming homoscedasticity}
#' }
