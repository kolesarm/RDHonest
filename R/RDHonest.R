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
#' @template RDse
#' @template RDoptBW
#' @template RDBW
#' @template RDclass
#' @template RDweights
#' @template Kern
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
#' \cite{Armstrong, Timothy B., and Michal Kolesár. 2018.
#' "Optimal Inference in a Class of Regression Models." Econometrica 86 (2):
#' 655–83.}
#'
#' \cite{Armstrong, Timothy B., and Michal Kolesár. 2020.
#' "Simple and Honest Confidence Intervals in Nonparametric Regression."
#' Quantitative Economics 11 (1): 1–39.}
#'
#' \cite{Imbens, Guido, and Kalyanaraman, Karthik,
#' "Optimal bandwidth choice for the regression discontinuity estimator." The
#' Review of Economic Studies 79 (3): 933-959.}
#'
#' \cite{Kolesár, Michal, and Christoph Rothe. 2018. "Inference in Regression
#' Discontinuity Designs with a Discrete Running Variable." American Economic
#' Review 108 (8): 2277–2304.}
#' }
#' @examples
#'
#' # Lee dataset
#' RDHonest(voteshare ~ margin, data = lee08, kern = "uniform", M = 0.1, h = 10)
#' RDHonest(cn~retired | elig_year, data=rcp, cutoff=0, M=c(4, 0.4),
#'           kern="triangular", opt.criterion="MSE", T0=0, h=20)
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
        message("Using ROT for smoothness constant, setting to ", M)
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
    ret$data <- d
    ret$call <- cl
    ret$na.action <- attr(mf, "na.action")

    ret
}
