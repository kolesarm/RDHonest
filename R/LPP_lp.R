#' Honest inference at a point
#'
#' Calculate estimators and one- and two-sided CIs based on local polynomial
#' estimator under second-order Taylor or Hölder smoothness class.
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
#' @template LPPseInitial
#' @return Returns an object of class \code{"LPPResults"}. The function
#'     \code{print} can be used to obtain and print a summary of the results. An
#'     object of class \code{"LPPResults"} is a list containing the following
#'     components
#'
#'     \describe{
#'   \item{\code{estimate}}{Point estimate. This estimate is MSE-optimal if
#'                   \code{opt.criterion="MSE"}}
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
#' }
#' @seealso \code{\link{LPPOptBW}}
#' @examples
#'
#' # Lee dataset
#' LPPHonest(voteshare ~ margin, data = lee08, subset = margin>0,
#'           kern = "uniform", M = 0.1, h = 10, sclass = "T")
#' @export
LPPHonest <- function(formula, data, subset, point=0, M, kern="triangular",
                      na.action, opt.criterion, h, se.method="nn", alpha=0.05,
                      beta=0.8, J=3, sclass="H", order=1, se.initial="EHW") {

    ## construct model frame
    cl <- mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

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
#' @template LPPseInitial
#' @return Returns an object of class \code{"LPPBW"}. The function \code{print}
#'     can be used to obtain and print a summary of the results. An object of
#'     class \code{"LPPBW"} is a list containing the following components:
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
#' @examples
#'
#' # Lee dataset
#' LPPOptBW(voteshare ~ margin, data = lee08, subset=margin>0,
#'          kern = "uniform", M = 0.1, opt.criterion = "MSE", sclass = "H")
#' @export
LPPOptBW <- function(formula, data, subset, point=0, M, kern="triangular",
                    na.action, opt.criterion, alpha=0.05, beta=0.8, sclass="H",
                    order=1, se.initial="EHW") {

    ## construct model frame
    cl <- mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    d <- LPPData(mf, point)

    ret <- LPPOptBW.fit(d, M, kern, opt.criterion, alpha, beta, sclass,
                       order, se.initial=se.initial)
    class(ret) <- "LPPBW"
    ret$call <- cl
    ret$na.action <- attr(mf, "na.action")

    ret
}

#' Honest inference at a point
#'
#' Basic computing engine called by \code{\link{LPPHonest}} to compute honest
#' confidence intervals for local polynomial estimators.
#' @param d object of class \code{"LPPData"}
#' @template RDse
#' @template RDoptBW
#' @template RDBW
#' @template RDclass
#' @template Kern
#' @template LPPseInitial
#' @return Returns an object of class \code{"LPPResults"}, see description in
#'     \code{\link{LPPHonest}}
## #' export
## LPPHonest.fit <- function(d, M, kern="triangular", h, opt.criterion,
##                          alpha=0.05, beta=0.8, se.method="nn",
##                          J=3, sclass="H", order=1, se.initial="EHW") {
##     CheckClass(d, "LPPData")

##     ## Initial se estimate
##     if (is.null(d$sigma2) & ("supplied.var" %in% se.method | missing(h)))
##         d <- NPRPrelimVar.fit(d, se.initial=se.initial)

##     if (missing(h)) {
##         r <- LPPOptBW.fit(d, M, kern, opt.criterion, alpha,
##                           beta, sclass, order)
##         h <- r$h
##     }

##     ## Suppress warnings about too few observations
##     r1 <- NPRreg.fit(d, h, kern, order, se.method, TRUE, J)
##     w <- r1$w(d$X)

##     ## If bandwidths too small
##     if (sum(w>0)==0) {
##         ## big bias / sd
##         bias <- sd <- upper <- hl <- sqrt(.Machine$double.xmax/10)
##         lower <- -upper
##     } else {
##         sd <- r1$se[se.method]
##         if(order==0) {
##             bias <- Inf
##         } else if (sclass=="T")  {
##             bias <- M/2 * sum(abs(w*d$X^2))
##             ## At boundary we know form of least favorable function
##         } else if (sclass=="H" & order==1 & length(unique(d$X>=0)==1)) {
##             bias <- abs(-M/2 * sum(w*d$X^2))
##         } else {
##             ww <- w[w>0]
##             xx <- d$X[w>0]
##             w2p <- function(s) abs(sum((ww*(xx-s))[xx>=s]))
##             w2m <- function(s) abs(sum((ww*(s-xx))[xx<=s]))
##             bp <- integrate(function(s) vapply(s, w2p, numeric(1)), 0, h)$value
##             bm <- integrate(function(s) vapply(s, w2m, numeric(1)), -h, 0)$value
##             bias <- M*(bp+bm)
##         }
##         lower <- r1$estimate - bias - stats::qnorm(1-alpha)*sd
##         upper <- r1$estimate + bias + stats::qnorm(1-alpha)*sd
##         hl <- CVb(bias/sd, alpha)$cv*sd
##     }

##     ## Finally, calculate coverage of naive CIs
##     z <- stats::qnorm(1-alpha/2)
##     naive <- stats::pnorm(z-bias/sd)-stats::pnorm(-z- bias/sd)

##     structure(list(estimate=r1$estimate, maxbias=bias, sd=sd,
##                    lower=lower, upper=upper, hl=hl, eff.obs=r1$eff.obs,
##                    h=h, naive=naive),
##               class="LPPResults")
## }


#' Optimal bandwidth selection for inference at a point
#'
#' Basic computing engine called by \code{\link{LPPOptBW}} used to find
#' optimal bandwidth
#' @param d object of class \code{"LPPData"}
#' @template RDoptBW
#' @template RDclass
#' @template Kern
#' @template LPPseInitial
#' @return a list with the following elements
#'     \describe{
#'     \item{\code{h}}{Bandwidth}
#'
#'     \item{\code{sigma2}}{estimate of conditional variance, from \code{d}}
#'    }
#' @examples
#' # Lee dataset
#' d <- LPPData(lee08[lee08$margin>0, ], point=0)
#' LPPOptBW.fit(d, kern = "uniform", M = 0.1, opt.criterion = "MSE")$h
#' @export
LPPOptBW.fit <- function(d, M, kern="triangular", opt.criterion, alpha=0.05,
                         beta=0.8, sclass="H", order=1, se.initial="EHW") {

    ## First check if sigma2 is supplied
    if (is.null(d$sigma2))
        d <- NPRPrelimVar.fit(d, se.initial=se.initial)

    ## Objective function for optimizing bandwidth
    obj <- function(h) {
        r <- NPRHonest.fit(d, M, kern, abs(h),
                          alpha=alpha, se.method="supplied.var",
                          sclass=sclass, order=order)
        if (opt.criterion=="OCI") {
            2*r$maxbias+r$sd*(stats::qnorm(1-alpha)+stats::qnorm(beta))
        } else if (opt.criterion=="MSE") {
            r$maxbias^2+r$sd^2
        } else if (opt.criterion=="FLCI") {
            r$hl
        } else {
            stop(sprintf("optimality criterion %s not yet implemented",
                         opt.criterion))
        }
    }

    hmin <- sort(unique(abs(d$X)))[order+1]
    hmax <- max(abs(d$X))
    ## Optimize piecewise constant function using modification of golden
    ## search. In fact, the criterion may not be unimodal, so proceed with
    ## caution. (For triangular kernel, it appears unimodal)
    if (kern=="uniform") {
        supp <- sort(unique(abs(d$X)))
        h <- gss(obj, supp[supp>=hmin])
    } else {
        h <- abs(stats::optimize(obj, interval=c(hmin, hmax),
                                 tol=.Machine$double.eps^0.75)$minimum)
    }

    list(h=h, sigma2=d$sigma2)
}


#' Rule of thumb bandwidth for inference at a point
#'
#' Calculate bandwidth for inference at a point on local linear regression using
#' method in Fan and Gijbels (1996, Chapter 4.2).
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
    X <- d$X

    if(is.null(boundary))
        boundary <- if ((min(X)>=0) | (max(X)<=0)) TRUE else FALSE
    if( (boundary==TRUE) & (order %% 2 ==0) )
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


#' @export
print.LPPBW <- function(x, digits = getOption("digits"), ...) {
    cat("Call:\n", deparse(x$call), "\n\n", sep = "", fill=TRUE)
    cat("Bandwidth: ", format(x$h, digits=digits), "\n\n")
    invisible(x)
}


#' @export
print.LPPResults <- function(x, digits = getOption("digits"), ...) {
    if (!is.null(x$call))
        cat("Call:\n", deparse(x$call), "\n\n", sep = "", fill=TRUE)

    cat("Inference by se.method:\n")
    r <- as.data.frame(x[c("estimate", "maxbias", "sd",
                           "lower", "upper", "hl")])
    names(r)[1:3] <- c("Estimate", "Maximum Bias", "Std. Error")
    r$name <- rownames(r)
    print.data.frame(as.data.frame(r[, 1:3]), digits=digits)
    cat("\nConfidence intervals:\n")
    fmt <- function(x) format(x, digits=digits, width=digits+1)

    for (j in seq_len(nrow(r)))
        cat(format(r[j, 7], width=5), " (", fmt(r[j, 1]-r[j, 6]),
            ", ", fmt(r[j, 1]+r[j, 6]), "), (", fmt(r[j, 4]),
            ", Inf), (-Inf, ", fmt(r[j, 5]), ")\n", sep="")

    cat("\nBandwidth: ", format(x$h, digits=digits), sep="")
    cat("\nNumber of effective observations:",
        format(x$eff.obs, digits=digits), "\n")

    invisible(x)
}
