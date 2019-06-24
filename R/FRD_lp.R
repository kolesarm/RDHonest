#' Honest inference in fuzzy RD
#'
#' Calculate estimators and one- and two-sided CIs based on local polynomial
#' estimator in fuzzy RD under second-order Taylor or Hölder smoothness class.
#'
#' The bandwidth is calculated to be optimal for a given performance criterion,
#' as specified by \code{opt.criterion}. It is calculated using the function
#' \code{\link{FRDOptBW}}. Alternatively, the bandwidth can be specified by
#' \code{h}.
#'
#' @template RDFormula
#' @template RDse
#' @template RDoptBW
#' @template RDBW
#' @template RDclass
#' @template Kern
#' @template bwequal
#' @template RDseInitial
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
#'   \item{\code{hp}, \code{hm}}{Bandwidths used above and below the cutoff}
#'
#'   \item{\code{naive}}{Coverage of CI that ignores bias and uses
#'                \code{qnorm(1-alpha/2)} as critical value}
#'
#'   \item{\code{call}}{the matched call}
#'
#'   \item{\code{fs}}{Estimate of the first-stage coefficient}
#'
#' }
#' @seealso \code{\link{FRDOptBW}}
#' @references{
#'
#' \cite{Armstrong, Tim, and Michal Kolesár. 2018. "Optimal Inference in a Class
#' of Regression Models." Econometrica 86 (2): 655–83.}
#' }
#' @examples
#' FRDHonest(cn~retired | elig_year, data=rcp, cutoff=0, M=c(1, 0.1),
#'           kern="triangular", opt.criterion="MSE", T0=0)
#' @export
FRDHonest <- function(formula, data, subset, cutoff=0, M, kern="triangular",
                     na.action, opt.criterion, bw.equal=TRUE, h,
                     se.method="nn", alpha=0.05, beta=0.8, J=3, sclass="H",
                     order=1, se.initial="EHW", T0=0) {

    ## construct model frame
    cl <- mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]

    formula <- Formula::as.Formula(formula)
    stopifnot(length(formula)[1] == 1L, length(formula)[2] == 2)
    mf$formula <- formula

    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    d <- FRDData(mf, cutoff)

    if (!missing(h)) {
        ret <- NPRHonest.fit(d, M, kern, h, alpha=alpha,
                            se.method=se.method, J=J, sclass=sclass,
                            order=order, se.initial=se.initial, T0=T0)
    } else {
        ret <- NPRHonest.fit(d, M, kern, opt.criterion=opt.criterion,
                            bw.equal=bw.equal, alpha=alpha, beta=beta,
                            se.method=se.method, J=J,
                            sclass=sclass, order=order, se.initial=se.initial, T0=T0)
    }

    ret$call <- cl
    ret$na.action <- attr(mf, "na.action")

    ret
}


#' Optimal Bandwidth Selection in Regression Discontinuity
#'
#' Estimate bandwidth for sharp RD based on local polynomial regression that
#' optimizes either maximum mean squared error, or length or quantiles of excess
#' length of a honest CI under second order Hölder or Taylor class.
#'
#' @template RDFormula
#' @template RDoptBW
#' @template RDclass
#' @template Kern
#' @template bwequal
#' @template RDseInitial
#' @param T0 Initial estimate of the treatment effect for calculating the
#'     optimal bandwidth. Only relevant for Fuzzy RD.
#' @return Returns an object of class \code{"RDBW"}. The function \code{print}
#'     can be used to obtain and print a summary of the results. An object of
#'     class \code{"RDBW"} is a list containing the following components:
#'
#'     \describe{
#'     \item{\code{hp}}{bandwidth for observations above cutoff}
#'
#'     \item{\code{hm}}{bandwidth for observations below cutoff, equal to \code{hp}
#'     unless \code{bw.equal==FALSE}}
#'
#'     \item{\code{sigma2m}, \code{sigma2p}}{estimate of conditional variance
#'     just above and just below cutoff, \eqn{\sigma^2_+(0)} and
#'     \eqn{\sigma^2_{-}(0)}}
#'
#'    \item{\code{f0}}{estimate of density of running variable at cutoff, if
#'    bandwidth computed using asymptotic method}
#'
#'    \item{\code{call}}{the matched call}
#'
#'    \item{\code{na.action}}{(where relevant) information on handling of missing
#'    data.}
#'
#'    }
#' @seealso \code{\link{RDHonest}}
#' @references{
#' \cite{Armstrong, Tim, and Michal Kolesár. 2018. "Optimal Inference in a Class
#' of Regression Models." Econometrica 86 (2): 655–83.}
#'
#' \cite{Armstrong, Timothy B., and Michal Kolesár. 2019.
#' "Simple and Honest Confidence Intervals in Nonparametric Regression", arXiv:
#' 1606.01200.}
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
#' FRDOptBW(cn~retired | elig_year, data=rcp, cutoff=0, M=c(1, 0.1),
#'           kern="triangular", opt.criterion="FLCI")
#' @export
FRDOptBW <- function(formula, data, subset, cutoff=0, M, kern="triangular",
                    na.action, opt.criterion, bw.equal=TRUE,
                    alpha=0.05, beta=0.8, sclass="H", order=1,
                    se.initial="EHW", T0=0) {

    ## construct model frame
    cl <- mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]

    formula <- Formula::as.Formula(formula)
    stopifnot(length(formula)[1] == 1L, length(formula)[2] == 2)
    mf$formula <- formula

    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    d <- FRDData(mf, cutoff)

    ret <- NPROptBW.fit(d, M, kern, opt.criterion, bw.equal, alpha, beta, sclass,
                       order, se.initial=se.initial, T0=T0)
    ret$call <- cl
    ret$na.action <- attr(mf, "na.action")

    ret
}
