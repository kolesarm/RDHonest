#' @param se.initial Method for estimating initial variance for computing
#'     optimal bandwidth. Except for \code{"nn"}, all methods assume
#'     homoskedasticity on either side of cutoff (for RD), or for all data (for
#'     inference at a point).
#'
#' \describe{
#'
#' \item{"EHW"}{Based on residuals from a local linear regression using a
#'             triangular kernel, and a bandwidth given by a rule-of-thumb
#'             bandwidth suggested by Fan and Gijbels (1996) (for inference at a
#'             point), or Imbens and Kalyanaraman (2012, IK) bandwidth (for
#'             fuzzy and sharp RD). For fuzzy RD, the IK bandwidth is based on
#'             the reduced-form regression.}
#'
#' \item{"Silverman"}{Use residuals from local constant regression with uniform
#' kernel and bandwidth selected using Silverman's rule of thumb, as in Equation
#' (14) in Imbens and Kalyanaraman (2012)}
#'
#' }
