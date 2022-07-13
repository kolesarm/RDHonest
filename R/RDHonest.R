#' Honest inference in RD
#'
#' Calculate estimators and bias-aware CIs for the sharp or fuzzy RD parameter,
#' or for value of the conditional mean at a point.
#'
#' The bandwidth is calculated to be optimal for a given performance criterion,
#' as specified by \code{opt.criterion}. Alternatively, for local polynomial
#' estimators, the bandwidth can be specified by \code{h}. For
#' \code{kern="optimal"}, calculate optimal estimators under second-order Taylor
#' smoothness class (sharp RD only).
#'
#' @template RDFormula
#' @param cutoff specifies the RD cutoff in the running variable. For inference
#'     at a point, specifies the point \eqn{x_0} at which to calculate the
#'     conditional mean.
#' @param kern specifies kernel function used in the local regression. It can
#'     either be a string equal to \code{"triangular"} (\eqn{k(u)=(1-|u|)_{+}}),
#'     \code{"epanechnikov"} (\eqn{k(u)=(3/4)(1-u^2)_{+}}), or \code{"uniform"}
#'     (\eqn{k(u)= (|u|<1)/2}), or else a kernel function. If equal to
#'     \code{"optimal"}, use the finite-sample optimal linear estimator under
#'     Taylor smoothness class, instead of a local linear estimator.
#' @param se.method Vector with methods for estimating standard error of
#'     estimate. If \code{NULL}, standard errors are not computed. The elements
#'     of the vector can consist of the following methods:
#'
#' \describe{
#'     \item{"nn"}{Nearest neighbor method}
#'
#'     \item{"EHW"}{Eicker-Huber-White, with residuals from local regression
#'     (local polynomial estimators only).}
#'
#'    \item{"supplied.var"}{Use conditional variance supplied by \code{sigma2}
#'         or \code{d} instead of computing residuals}
#'
#' }
#' @param J Number of nearest neighbors, if "nn" is specified in
#'     \code{se.method}.
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
#' @param M Bound on second derivative of the conditional mean function.
#' @param sclass Smoothness class, either \code{"T"} for Taylor or
#'     \code{"H"} for Hölder class.
#' @param h bandwidth, a scalar parameter. If not supplied, optimal bandwidth is
#'     computed according to criterion given by \code{opt.criterion}.
#' @param weights Optional vector of weights to weight the observations
#'     (useful for aggregated data). Disregarded if optimal kernel is used.
#' @param point.inference Do inference at a point determined by \code{cutoff}
#'     instead of RD.
#' @param T0 Initial estimate of the treatment effect for calculating the
#'     optimal bandwidth. Only relevant for Fuzzy RD.
#' @return Returns an object of class \code{"NPRResults"}. The function
#'     \code{print} can be used to obtain and print a summary of the results. An
#'     object of class \code{"NPRResults"} is a list containing the following
#'     components
#'
#'     \describe{
#'   \item{\code{estimate}}{Point estimate. This estimate is MSE-optimal if
#'                   \code{opt.criterion="MSE"}}
#'
#'   \item{\code{lff}}{Least favorable function, only relevant for optimal
#'              estimator under Taylor class.}
#'
#'   \item{\code{maxbias}}{Maximum bias of \code{estimate}}
#'
#'   \item{\code{sd}}{Standard deviation of estimate}
#'
#'   \item{\code{lower}, \code{upper}}{Lower (upper) end-point of a one-sided CI
#'         based on \code{estimate}. This CI is optimal if
#'         \code{opt.criterion=="OCI"}}
#'
#'   \item{\code{hl}}{Half-length of a two-sided CI based on \code{estimate}, so
#'             that the CI is given by \code{c(estimate-hl, estimate+hl)}. The
#'             CI is optimal if \code{opt.criterion="FLCI"}}
#'
#'   \item{\code{eff.obs}}{Effective number of observations used by
#'             \code{estimate}}
#'
#'   \item{\code{h}}{Bandwidth used}
#'
#'   \item{\code{naive}}{Coverage of CI that ignores bias and uses
#'                \code{qnorm(1-alpha/2)} as critical value}
#'
#'   \item{\code{call}}{the matched call}
#'
#'   \item{\code{fs}}{Estimate of the first-stage coefficient (sharp RD only)}
#'
#' }
#' @references{
#'
#' \cite{Timothy B. Armstrong and Michal Kolesár. Optimal inference in a class
#' of regression models. Econometrica, 86(2):655–683, March 2018.
#' \doi{10.3982/ECTA14434}}
#'
#' \cite{Timothy B. Armstrong and Michal Kolesár. Simple and honest confidence
#' intervals in nonparametric regression. Quantitative Economics, 11(1):1–39,
#' January 2020. \doi{10.3982/QE1199}}
#'
#' \cite{Guido W. Imbens and Karthik Kalyanaraman. Optimal bandwidth choice for
#' the regression discontinuity estimator. The Review of Economic Studies,
#' 79(3):933–959, July 2012. \doi{10.1093/restud/rdr043}}
#'
#' \cite{Michal Kolesár and Christoph Rothe. Inference in regression
#' discontinuity designs with a discrete running variable. American Economic
#' Review, 108(8):2277—-2304, August 2018. \doi{10.1257/aer.20160945}}
#'
#' }
#' @examples
#'
#' # Lee dataset
#' RDHonest(voteshare ~ margin, data = lee08, kern = "uniform", M = 0.1, h = 10)
#' RDHonest(cn~retired | elig_year, data=rcp, cutoff=0, M=c(4, 0.4),
#'           kern="triangular", opt.criterion="MSE", T0=0, h=3)
#' RDHonest(voteshare ~ margin, data = lee08, subset = margin>0,
#'           kern = "uniform", M = 0.1, h = 10, point.inference=TRUE)
#' @export
RDHonest <- function(formula, data, subset, weights, cutoff=0, M,
                     kern="triangular", na.action, opt.criterion="MSE", h,
                     se.method="nn", alpha=0.05, beta=0.8, J=3, sclass="H",
                     T0=0, point.inference=FALSE) {

    ## construct model frame
    cl <- mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    formula <- Formula::as.Formula(formula)
    ## one LHS, at most 2 RHS
    stopifnot(length(formula)[1] == 1L, length(formula)[2] <= 2)
    mf$formula <- formula

    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mf$weights  <- mf$"(weights)"
    d <- if (point.inference) {
             LPPData(mf, cutoff)
         } else if (length(formula)[2]==2) {
             FRDData(mf, cutoff)
         } else {
             RDData(mf, cutoff)
         }

    if (missing(M)) {
        M <- NPR_MROT.fit(d)
        message("Using Armstong Kolesar (2020) ROT for smoothness constant M")
    }
    if (kern=="optimal") {
        ret <- RDTOpt.fit(d, M, opt.criterion, alpha, beta, se.method, J)
    } else if (!missing(h)) {
        ret <- NPRHonest.fit(d, M, kern, h, alpha=alpha, se.method=se.method,
                             J=J, sclass=sclass, T0=T0)
    } else {
        ret <- NPRHonest.fit(d, M, kern, opt.criterion=opt.criterion,
                             alpha=alpha, beta=beta, se.method=se.method, J=J,
                             sclass=sclass, T0=T0)
    }
    ret$call <- cl
    ret$na.action <- attr(mf, "na.action")

    ret
}
