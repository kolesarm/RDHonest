#' Honest inference at a point
#'
#' Calculate estimators and one- and two-sided honest CIs for value of
#' conditional mean at a point based on a local polynomial estimator under
#' second-order Taylor or Hölder smoothness class.
#'
#' The bandwidth is calculated to be optimal for a given performance criterion,
#' as specified by \code{opt.criterion}. It is calculated using the function
#' \code{\link{LPPOptBW}}. Alternatively, the bandwidth can be specified by
#' \code{h}.
#'
#' @template LPPFormula
#' @template RDse
#' @template RDoptBW
#' @template RDBW
#' @template RDclass
#' @template Kern
#' @template RDseInitial
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
#' @seealso \code{\link{LPPOptBW}}
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
                      se.method="nn", alpha=0.05, beta=0.8, J=3, sclass="H",
                      order=1, se.initial="EHW") {

    ## construct model frame
    cl <- mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mf$weights  <- mf$"(weights)"
    d <- LPPData(mf, point)

    if (!missing(h)) {
        ret <- NPRHonest.fit(d, M, kern, h, alpha=alpha,
                             se.method=se.method, J=J, sclass=sclass,
                             order=order, se.initial=se.initial)
    } else {
        ret <- NPRHonest.fit(d, M, kern, opt.criterion=opt.criterion,
                             alpha=alpha, beta=beta, se.method=se.method, J=J,
                             sclass=sclass, order=order, se.initial=se.initial)
    }

    ret$call <- cl
    ret$na.action <- attr(mf, "na.action")

    ret
}


#' Optimal Bandwidth Selection for inference at a point
#'
#' Estimate bandwidth based on local polynomial regression that optimizes either
#' maximum mean squared error, or length or quantiles of excess length of a
#' honest CI under second order Hölder or Taylor class.
#'
#' @template LPPFormula
#' @template RDoptBW
#' @template RDclass
#' @template Kern
#' @template RDseInitial
#' @return Returns an object of class \code{"NPRBW"}. The function \code{print}
#'     can be used to obtain and print a summary of the results. An object of
#'     class \code{"NPRBW"} is a list containing the following components:
#'
#'     \describe{
#'     \item{\code{h}}{Bandwidth}
#'
#'     \item{\code{sigma2}}{estimate of conditional variance at a point}
#'
#'    \item{\code{call}}{the matched call}
#'
#'    \item{\code{na.action}}{(where relevant) information on handling of
#'                            missing data.}
#'
#'    }
#' @seealso \code{\link{LPPHonest}}
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
#' LPPOptBW(voteshare ~ margin, data = lee08, subset=margin>0,
#'          kern = "uniform", M = 0.1, opt.criterion = "MSE", sclass = "H")
#' @export
LPPOptBW <- function(formula, data, subset, weights, point=0, M,
                     kern="triangular", na.action, opt.criterion, alpha=0.05,
                     beta=0.8, sclass="H", order=1, se.initial="EHW") {

    ## construct model frame
    cl <- mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mf$weights  <- mf$"(weights)"
    d <- LPPData(mf, point)

    ret <- NPROptBW.fit(d, M, kern, opt.criterion, alpha, beta, sclass,
                       order, se.initial=se.initial)
    ret$call <- cl
    ret$na.action <- attr(mf, "na.action")

    ret
}


#' Rule of thumb bandwidth for inference at a point
#'
#' Calculate bandwidth for inference at a point with local linear regression
#' using method in Fan and Gijbels (1996, Chapter 4.2).
#'
#' @param d object of class \code{"LPPData"}
#' @template Kern
#' @return ROT bandwidth
#' @param boundary Is point at a boundary?
#' @references{
#'
#' \cite{Fan , J., and I. Gijbels (1996): Local Polynomial Modelling and Its
#' Applications, Monographs on Statistics and Applied Probability. Chapman &
#' Hall/CRC, New York, NY.}
#'
#' }
#' @examples
#' dp <- LPPData(lee08[lee08$margin>0, ], point=0)
#' bp1 <- ROTBW.fit(dp, kern="uniform", order=1)
#' @export
ROTBW.fit <- function(d, kern="triangular", order=1, boundary=NULL) {
    CheckClass(d, "LPPData")
    X <- d$X

    if(is.null(boundary))
        boundary <- if ((min(X)>=0) | (max(X)<=0)) TRUE else FALSE
    if((boundary==TRUE) & (order %% 2 ==0))
        warning("ROT method for computing bandwidth requires either\n",
                "order to be odd or else a boundary point")
    N <- length(d$X)

    ## STEP 0: Estimate f_X(0) using Silverman
    h1 <- 1.843 *
        min(stats::sd(X), (stats::quantile(X, 0.75) -
                           stats::quantile(X, 0.25)) / 1.349) / N^(1/5)
    f0 <- sum(abs(X) <= h1) / (2*N*h1)

    ## STEP 1: Estimate (p+1)th derivative and sigma^2 using global polynomial
    ## regression
    r1 <- stats::lm(d$Y ~ 0 + outer(X, 0:(order+3), "^"))
    deriv <- unname(r1$coefficients[order+2])
    sigma2 <- stats::sigma(r1)^2

    ## STEP 2: Kernel constants
    if (is.function(kern)) {
        ke <- EqKern(kern, boundary=boundary, order=order)
        nu0 <- KernMoment(ke, moment=0, boundary=boundary, "raw2")
        mup <- KernMoment(ke, moment=order+1, boundary=boundary, "raw")
    } else {
        s <- RDHonest::kernC[RDHonest::kernC$kernel==kern &
                             RDHonest::kernC$order==order &
                             RDHonest::kernC$boundary==boundary, ]
        nu0 <- s$nu0
        mup <- s[[paste0("mu", order+1)]]
    }

    ## STEP 3: Plug in
    B <- deriv * mup
    V <- sigma2 * nu0 /f0

    (V/(B^2 * 2 * (order+1) * N))^(1/(2*order+3))
}
