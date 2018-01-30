#' Local Polynomial Regression at a point
#'
#' Calculate estimate and its variance given a bandwidth using local polynomial
#' regression of order \code{order}.
#'
#' @param d object of class \code{"LPPData"}
#' @template LPPBW
#' @template Kern
#' @template RDse
#' @param no.warning Don't warn about too few observations
#' @return list with elements:
#'
#' \describe{
#'     \item{estimate}{point estimate}
#'     \item{se}{Named vector of standard error estimates, as specified
#'               by \code{se.method}.}
#'     \item{w}{Implicit weight function used}
#'
#'     \item{sigma2}{Estimate of \eqn{sigma^2(X)} for values of \eqn{X}
#'             receiving positive kernel weight. By default, estimates are based
#'             on squared regression residuals, as used in \code{"EHW"}. If
#'             \code{"demeaned"} or \code{"nn"} is specifed, estimates are based
#'             on that method, with \code{"nn"} method used if both are
#'             specified.}
#'
#'      \item{eff.obs}{Number of effective observations}
#'
#' }
#' @export
LPPreg <- function(d, h, kern="triangular", order=1, se.method="nn",
                   no.warning=FALSE, J=3) {

    K <- if (!is.function(kern)) {
             EqKern(kern, boundary=FALSE, order=0)
         } else {
             kern
         }
    W <- if (h<=0) 0*d$X else K(d$X/h) # kernel weights

    ## variance calculations are faster if we only keep data with positive
    ## kernel weights
    d$X <- d$X[W>0]

    if ((length(d$X) < 3*order) & no.warning==FALSE) {
        warning("Too few observations to compute estimates.\nOnly ",
                length(d$X), " units with positive weights")
    }
    if (length(unique(d$X)) <= order & no.warning==FALSE) {
        warning("Too few distinct values to compute estimates.\nOnly ",
                length(unique(d$X)), " unique values for ",
                "independent variable with positive weights")
    }

    r <- LPReg(d$X, d$Y[W>0], h, K, order, se.method, d$sigma2[W>0], J)

    ## Plug-in Asymptotic variance not computed
    plugin <- NA
    names(plugin) <- "plugin"
    se <- sqrt(c(r$var, plugin))
    list(estimate=r$theta, se=se, w=r$w, sigma2=r$sigma2, eff.obs=r$eff.obs)
}


#' Compute preliminary estimate of variance
#'
#' Compute estimate of variance, which can then be used in optimal bandwidth
#' calculations.
#'
#' @param d object of class \code{"LPPData"}
#' @template LPPseInitial
#' @return object of class \code{"LPPData"} containing estimated variances.
#' @export
LPPPrelimVar <- function(d, se.initial="IKEHW") {
    if (se.initial=="ROTdemeaned") {
        r1 <- LPPreg(d, ROTBW.fit(d), se.method="demeaned")
        d$sigma2 <- rep(mean(r1$sigma2), length(d$X))
    } else if (se.initial=="ROTEHW"){
        r1 <- LPPreg(d, ROTBW.fit(d), se.method="EHW")
        d$sigma2 <- rep(mean(r1$sigma2), length(d$X))
    } else {
        stop("Unknown method for estimating initial variance")
    }

    d
}
