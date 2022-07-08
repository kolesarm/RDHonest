#' Honest inference at a point
#'
#' Calculate estimators and one- and two-sided honest CIs for value of
#' conditional mean at a point based on a local polynomial estimator under
#' second-order Taylor or Hölder smoothness class.
#'
#' The bandwidth is calculated to be optimal for a given performance criterion,
#' as specified by \code{opt.criterion}. Alternatively, the bandwidth can be
#' specified by \code{h}.
#'
#' @template LPPFormula
#' @template RDse
#' @template RDoptBW
#' @template RDBW
#' @template RDclass
#' @template Kern
#' @return Returns an object of class \code{"NPResults"}. The function
#'     \code{print} can be used to obtain and print a summary of the results. An
#'     object of class \code{"NPRResults"} is a list containing the following
#'     components
#'
#'     \describe{
#'   \item{\code{estimate}}{Point estimate. This estimate is MSE-optimal if
#'                   \code{opt.criterion="MSE"}}
#'
#'   \item{\code{lff}}{Not relevant for inference at a point}
#'
#'   \item{\code{maxbias}}{Maximum bias of \code{estimate}}
#'
#'   \item{\code{sd}}{Standard deviation of estimate}
#'
#'   \item{\code{lower}, \code{upper}}{Lower (upper) end-point of a one-sided CI
#'         based on \code{estimate}. This CI is optimal if
#'         \code{opt.criterion="OCI"}}
#'
#'   \item{\code{hl}}{Half-length of a two-sided CI based on \code{estimate}, so
#'          the CI is \code{c(estimate-hl, estimate+hl)}. The CI is optimal if
#'          \code{opt.criterion="FLCI"}}
#'
#'   \item{\code{eff.obs}}{Effective number of observations used by
#'             \code{estimate}}
#'
#'   \item{\code{h}}{Bandwidth used}
#'
#'   \item{\code{naive}}{Coverage of CI that ignores bias and uses
#'                \code{qnorm(1-alpha/2)} as critical value}
#'
#'   \item{\code{call}}{The matched call}
#'
#'   \item{\code{fs}}{Not relevant for inference at a point}
#'
#' }
#' @references{
#'
#' \cite{Armstrong, Timothy B., and Michal Kolesár. 2020.
#' "Simple and Honest Confidence Intervals in Nonparametric Regression."
#' Quantitative Economics 11 (1): 1–39.}
#'
#' }
#' @examples
#'
#' # Lee dataset
#' LPPHonest(voteshare ~ margin, data = lee08, subset = margin>0,
#'           kern = "uniform", M = 0.1, h = 10, sclass = "T")
#' @export
LPPHonest <- function(formula, data, subset, weights, point=0, M,
                      kern="triangular", na.action, opt.criterion, h,
                      se.method="nn", alpha=0.05, beta=0.8, J=3, sclass="H") {

    ## construct model frame
    cl <- mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)

    mf <- eval(mf, parent.frame())
    mf$weights  <- mf$"(weights)"
    d <- LPPData(mf, point)
    if (missing(M))
        M <- NPR_MROT.fit(d)

    if (!missing(h)) {
        ret <- NPRHonest.fit(d, M, kern, h, alpha=alpha,
                             se.method=se.method, J=J, sclass=sclass)
    } else {
        ret <- NPRHonest.fit(d, M, kern, opt.criterion=opt.criterion,
                             alpha=alpha, beta=beta, se.method=se.method, J=J,
                             sclass=sclass)
    }

    ret$call <- cl
    ret$na.action <- attr(mf, "na.action")

    ret
}
