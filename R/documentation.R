#' Lee (2008) US House elections dataset
#'
#' @format A data frame with 6558 rows and 2 variables:
#' \describe{
#'   \item{voteshare}{Vote share in next election}
#'   \item{margin}{Democratic margin of victory}
#'  }
#' @source Mostly Harmless Econometrics website
"lee08"

#' Constants for common kernels.
#'
#' First four moments of uniform, triangular, and Epanechnikov equivalent
#' kernels. Up to numerical integration precision, these moments are matched by
#' \code{KernMoment()}
#'
#' @format A data frame with 18 rows and 19 variables:
#' \describe{
#'   \item{kernel}{Kernel type.}
#'   \item{order}{Order of local polynomial.}
#'   \item{boundary}{Boundary regression?}
#'   \item{mu0, mu1, mu2, mu3, mu4}{\eqn{\int_X u^j k(u)}, raw moments}
#'   \item{nu0, nu1, nu2, nu3, nu4}{\eqn{\int_X u^j k^2(u)}, raw moments of
#'         kernel squared}
#'   \item{pi0, pi1, pi2, pi3, pi4}{\eqn{\int_X abs{u^j k(u)}}, absolute moments}
#'   \item{pMSE}{constant for pointwise MSE optimal bandwidth,
#'        \eqn{((p+1)!^2\nu_0 / (2(p+1)\mu_{p+1}^2))^{1/(2p+3)}}, see page 67 in
#'        Fan and Gijbels}
#'  }
#' @source Computed analytically using symbolic math software
"kernC"
