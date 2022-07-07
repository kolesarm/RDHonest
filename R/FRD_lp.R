#' Honest inference in fuzzy RD
#'
#' Calculate estimators and one- and two-sided CIs based on local polynomial
#' estimator in fuzzy RD under second-order Taylor or Hölder smoothness class.
#'
#' The bandwidth is calculated to be optimal for a given performance criterion,
#' as specified by \code{opt.criterion}. Alternatively, the bandwidth can be
#' specified by \code{h}.
#'
#' @template RDFormula
#' @template RDse
#' @template RDoptBW
#' @template RDBW
#' @template RDweights
#' @template RDclass
#' @template Kern
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
#'   \item{\code{lff}}{Not relevant for fuzzy RD.}
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
#'   \item{\code{fs}}{Estimate of the first-stage coefficient}
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
#' }
#' @examples
#' FRDHonest(cn~retired | elig_year, data=rcp, cutoff=0, M=c(5, 0.5),
#'           kern="triangular", opt.criterion="MSE", T0=0)
#' @export
FRDHonest <- function(formula, data, subset, weights, cutoff=0, M,
                      kern="triangular", na.action, opt.criterion, h,
                      se.method="nn", alpha=0.05, beta=0.8, J=3, sclass="H",
                      order=1, T0=0) {

    ## construct model frame
    cl <- mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]

    formula <- Formula::as.Formula(formula)
    stopifnot(length(formula)[1] == 1L, length(formula)[2] == 2)
    mf$formula <- formula

    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mf$weights  <- mf$"(weights)"
    d <- FRDData(mf, cutoff)
    if (missing(M))
        M <- NPR_MROT.fit(d)

    if (!missing(h)) {
        ret <- NPRHonest.fit(d, M, kern, h, alpha=alpha, se.method=se.method,
                             J=J, sclass=sclass, order=order, T0=T0)
    } else {
        ret <- NPRHonest.fit(d, M, kern, opt.criterion=opt.criterion,
                             alpha=alpha, beta=beta, se.method=se.method, J=J,
                             sclass=sclass, order=order, T0=T0)
    }

    ret$call <- cl
    ret$na.action <- attr(mf, "na.action")

    ret
}
