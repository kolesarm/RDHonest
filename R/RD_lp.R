#' Honest inference in RD
#'
#' Calculate estimators and bias-aware one- and two-sided CIs for the sharp RD
#' parameter.
#'
#' The bandwidth is calculated to be optimal for a given performance criterion,
#' as specified by \code{opt.criterion}. Alternatively, for local polynomial
#' estimators, the bandwidth can be specified by \code{h}. If
#' \code{kern="optimal"}, calculate optimal estimators under second-order Taylor
#' smoothness class.
#'
#' @template RDFormula
#' @template RDse
#' @template RDoptBW
#' @template RDBW
#' @template RDclass
#' @template Kern
#' @template RDseInitial
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
#'   \item{\code{fs}}{Not relevant for sharp RD}
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
#' @export
RDHonest <- function(formula, data, subset, weights, cutoff=0, M,
                     kern="triangular", na.action, opt.criterion="MSE", h,
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
    d <- RDData(mf, cutoff)
    if (missing(M))
        M <- NPR_MROT.fit(d)

    if (kern=="optimal") {
        ret <- RDTOpt.fit(d, M, opt.criterion=opt.criterion, alpha=alpha,
                          beta=beta, se.method=se.method, J=J,
                          se.initial=se.initial)
    } else if (!missing(h)) {
        ret <- NPRHonest.fit(d, M, kern, h, alpha=alpha, se.method=se.method,
                             J=J, sclass=sclass, order=order,
                             se.initial=se.initial)
    } else {
        ret <- NPRHonest.fit(d, M, kern, opt.criterion=opt.criterion,
                             alpha=alpha, beta=beta, se.method=se.method, J=J,
                             sclass=sclass, order=order, se.initial=se.initial)
    }

    ret$call <- cl
    ret
}


## Imbens and Kalyanaraman bandwidth. Only used by NPRPrelimVar.fit
##
##  Reproduce bandwidth from Section 6.2 in Imbens and Kalyanaraman (2012)
IKBW.fit <- function(d, kern="triangular", order=1, verbose=FALSE) {
    if (order!=1)
        stop("Only works for local linear regression.")
    CheckClass(d, "RDData")

    X <- c(d$Xm, d$Xp)
    Nm <- length(d$Xm)
    Np <- length(d$Xp)
    N <- Nm+Np

    ## STEP 0: Kernel constant
    if (is.character(kern)) {
        s <- RDHonest::kernC[with(RDHonest::kernC,
                        order==1 & boundary==TRUE & kernel==kern, ), ]
    } else if (is.function(kern)) {
        ke <- EqKern(kern, boundary=TRUE, order=1)
        s <- list(mu2=KernMoment(ke, moment=2, boundary=TRUE, "raw"),
                  nu0=KernMoment(ke, moment=0, boundary=TRUE, "raw2"))
    }
    const <- (s$nu0/s$mu2^2)^(1/5)

    ## STEP 1: Estimate f(0), sigma^2_(0) and sigma^2_+(0), using Silverman
    ## pilot bandwidth for uniform kernel
    d <- NPRPrelimVar.fit(d, se.initial="Silverman")
    h1 <- 1.84*stats::sd(X)/N^(1/5)
    f0 <- sum(abs(X) <= h1) / (2*N*h1)
    varm <- d$sigma2m[1]
    varp <- d$sigma2p[1]

    ## STEP 2: Estimate second derivatives m_{+}^(2) and m_{-}^(2)

    ## Estimate third derivative using 3rd order polynomial: Equation (14)
    m3 <- 6*stats::coef(stats::lm(I(c(d$Ym, d$Yp)) ~ I(X>=0) + X +
                                      I(X^2) + I(X^3)))[5]

    ## Left and right bandwidths, Equation (15) and page 956.
    ## Optimal constant based on one-sided uniform Kernel, 7200^(1/7),
    h2m <- 7200^(1/7) * (varm/(f0*m3^2))^(1/7) * Nm^(-1/7)
    h2p <- 7200^(1/7) * (varp/(f0*m3^2))^(1/7) * Np^(-1/7)

    ## estimate second derivatives by local quadratic
    m2m <- 2*stats::coef(stats::lm(d$Ym ~ d$Xm + I(d$Xm^2),
                                   subset=(d$Xm >= -h2m)))[3]
    m2p <- 2*stats::coef(stats::lm(d$Yp ~ d$Xp + I(d$Xp^2),
                                   subset=(d$Xp <= h2p)))[3]

    ## STEP 3: Calculate regularization terms, Equation (16)
    rm <- 2160*varm / (sum(d$Xm >= -h2m) * h2m^4)
    rp <- 2160*varp / (sum(d$Xp <= h2p) * h2p^4)

    if(verbose)
        cat("\n h1: ", h1, "\n N_{-}, N_{+}: ", Nm, Np, "\n f(0): ", f0,
            "\n sigma^2_{+}(0): ", sqrt(varp),
            "^2\n sigma^2_{+}(0):", sqrt(varm), "^2",
            "\n m3: ", m3, "\n h_{2, +}:", h2p, "h_{2, -}:", h2m,
            "\n m^{(2)}_{+}: ", m2p, "m^{(2)}_{-}: ", m2m,
            "\n r_{+}:", rp, "\n r_{-}:", rm, "\n\n")

    ## Final bandwidth: Equation (17)
    unname(const * ((varp+varm) / (f0*N * ((m2p-m2m)^2+rm+rp)))^(1/5))
}
