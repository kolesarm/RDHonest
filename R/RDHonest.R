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
#'    \item{"supplied.var"}{Use conditional variance supplied by \code{sigmaY2}
#'         instead of computing residuals}
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
#'     The methods use conditional variance given by \code{sigmaY2}, if
#'     supplied. Otherwise, for the purpose of estimating the optimal bandwidth,
#'     conditional variance is estimated using the method specified by
#'     \code{se.initial}.
#' @param beta Determines quantile of excess length to optimize, if bandwidth
#'     optimizes given quantile of excess length of one-sided confidence
#'     intervals; otherwise ignored.
#' @param alpha determines confidence level, \code{1-alpha} for
#'     constructing/optimizing confidence intervals.
#' @param M Bound on second derivative of the conditional mean function.
#' @param sclass Smoothness class, either \code{"T"} for Taylor or \code{"H"}
#'     for Hölder class.
#' @param h bandwidth, a scalar parameter. If not supplied, optimal bandwidth is
#'     computed according to criterion given by \code{opt.criterion}.
#' @param weights Optional vector of weights to weight the observations (useful
#'     for aggregated data). Disregarded if optimal kernel is used.
#' @param point.inference Do inference at a point determined by \code{cutoff}
#'     instead of RD.
#' @param T0 Initial estimate of the treatment effect for calculating the
#'     optimal bandwidth. Only relevant for Fuzzy RD.
#' @param sigmaY2 Supply variance of outcome. Ignored when kernel is optimal.
#' @param sigmaD2 Supply variance of treatment (fuzzy RD only).
#' @param sigmaYD Supply covariance of treatment and outcome (fuzzy RD only).
#' @param clusterid Cluster id for cluster-robust standard errors
#' @return Returns an object of class \code{"RDResults"}. The function
#'     \code{print} can be used to obtain and print a summary of the results. An
#'     object of class \code{"RDResults"} is a list containing the following
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
#' RDHonest(voteshare ~ margin, data = lee08, kern = "uniform", M = 0.1, h = 10)
#' RDHonest(cn | retired ~ elig_year, data=rcp, cutoff=0, M=c(4, 0.4),
#'           kern="triangular", opt.criterion="MSE", T0=0, h=3)
#' RDHonest(voteshare ~ margin, data = lee08, subset = margin>0,
#'           kern = "uniform", M = 0.1, h = 10, point.inference=TRUE)
#' @export
RDHonest <- function(formula, data, subset, weights, cutoff=0, M,
                     kern="triangular", na.action, opt.criterion="MSE", h,
                     se.method="nn", alpha=0.05, beta=0.8, J=3, sclass="H",
                     T0=0, point.inference=FALSE, sigmaY2, sigmaD2, sigmaYD,
                     clusterid) {
    ## construct model frame
    cl <- mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", "sigmaY2",
                 "sigmaD2", "sigmaYD", "clusterid"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    formula <- Formula::as.Formula(formula)
    ## at most 2 LHS and at most 2 RHS
    stopifnot(length(formula)[1] <= 2L, length(formula)[2] <= 2L)

    if (point.inference) {
        method <- c("Value of conditional mean"="IP")
    } else if (length(formula)[1]==2) {
        method <- c("Fuzzy RD parameter"="FRD")
    } else {
        method <- c("Sharp RD parameter"="SRD")
    }
    mf$formula <- formula
    ## http://madrury.github.io/jekyll/update/statistics/2016/07/20/lm-in-R.html
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    d <- NPRData(mf, cutoff, method, formula)
    process_options(M, se.method, method, d, kern)

    if (!is.null(d$covs) && (missing(M) || missing(h))) {
        ## remove covariates and compute bw w/o covariates
        d0 <- d
        d0$covs <- NULL
        ret <- NPRHonest(d0, MROT(d0), kern, h, opt.criterion=opt.criterion,
                         alpha=alpha, beta=beta, se.method=se.method, J=J,
                         sclass=sclass, T0=T0)
        ## Compute covariate-adjusted outcome
        d <- covariate_adjust(d, kern, ret$coefficients$bandwidth)
    }

    if (missing(M)) {
        M <- MROT(d)
        message("Using Armstong & Kolesar (2020) ROT for smoothness constant M")
    }

    if (kernel_type(kern)=="optimal") {
        ret <- RDTOpt(d, M, opt.criterion, alpha, beta, se.method, J)
    } else {
        ret <- NPRHonest(d, M, kern, h, opt.criterion=opt.criterion,
                         alpha=alpha, beta=beta, se.method=se.method, J=J,
                         sclass=sclass, T0=T0)
    }

    ret$call <- ret$lm$call <- cl
    ret$na.action <- attr(mf, "na.action")

    if (!is.finite(ret$coefficients$leverage) || ret$coefficients$leverage>0.1)
        message(paste0("Maximal leverage is large: ",
                       round(ret$coefficients$leverage, 2),
                       ".\nInference may be inaccurate. ",
                       "Consider using bigger bandwidth."))
    ret$coefficients$term <- names(method)
    ret
}

covariate_adjust <- function(d, kern, h) {
    if (is.null(d$Y_unadj)) d$Y_unadj <- d$Y
    d0 <- d
    d0$Y <- d0$Y_unadj
    r <- NPReg(d0, h, kern, order=1, se.method="EHW")
    be <- as.matrix(r$lm$coefficients)
    L <- NCOL(d$covs)
    d$Y <-d$Y_unadj-d$covs %*% be[seq_len(L)+NROW(be)-L, , drop=FALSE]
    d
}



## @param T0 Initial estimate of the treatment effect for calculating the
##     optimal bandwidth. Only relevant for Fuzzy RD.
## @param T0bias When evaluating the maximum bias of the estimate, use the
##     estimate itself (if \code{T0bias==FALSE}), or use the preliminary
##     estimate \code{T0} (if \code{T0bias==TRUE}). Only relevant for Fuzzy RD.
NPRHonest <- function(d, M, kern="triangular", h, opt.criterion, alpha=0.05,
                      beta=0.8, se.method="nn", J=3, sclass="H", T0=0,
                      T0bias=FALSE) {
    if (missing(h)) {
        d0 <- d
        d0$covs <- d0$Y_unadj <- NULL
        h <- OptBW(d0, M, kern, opt.criterion, alpha, beta, sclass, T0)
    }
    if (!is.null(d$Y_unadj)) d$Y <- d$Y_unadj

    r1 <- NPReg(d, h, kern, order=1, se.method, J)
    wt <- r1$est_w[r1$est_w!=0]
    xx <- d$X[r1$est_w!=0]
    ## Are we at a boundary?
    bd <- TRUE
    if (inherits(d, "IP")) {
        bd <- length(unique(xx>=0))==1
    }

    if (T0bias && inherits(d, "FRD")) {
        ## multiply bias and sd by r1$fs to make if free of first stage
        r1$se <- r1$se*abs(r1$fs)
        M <- unname(c((M[1]+M[2]*abs(T0)), M))
    } else if (!T0bias && inherits(d, "FRD")) {
        M <- unname(c((M[1]+M[2]*abs(r1$estimate)) / abs(r1$fs), M))
    } else {
        M[2:3] <- c(NA, NA)
    }

    ## Determine bias
    if (r1$eff.obs==0) {
        ## If bandwidths too small, big bias / sd
        bias <- r1$se <- sqrt(.Machine$double.xmax/10)
    } else if (sclass=="T")  {
        bias <- M[1]/2 * (sum(abs(wt*xx^2)))
    } else if (sclass=="H" && bd) {
        ## At boundary we know form of least favorable function
        bias <- M[1]/2 *
            abs(sum(wt[xx<0]*xx[xx<0]^2)-sum(wt[xx>=0]*xx[xx>=0]^2))
    } else {
        ## Else need to find numerically (same formula if order=2)
        w2p <- function(s) abs(sum((wt * (xx-s))[xx>=s]))
        w2m <- function(s) abs(sum((wt * (s-xx))[xx<=s]))
        bp <- stats::integrate(function(s) vapply(s, w2p, numeric(1)), 0, h)
        bm <- stats::integrate(function(s) vapply(s, w2m, numeric(1)), -h, 0)
        bias <- M[1] * (bp$value+bm$value)
    }
    method <- switch(sclass, H="Holder", "Taylor")

    d$est_w <- r1$est_w
    d$sigma2 <- r1$sigma2
    co <- data.frame(term=NA, estimate=r1$estimate, std.error=r1$se,
                     maximum.bias=bias, conf.low=NA, conf.high=NA,
                     conf.low.onesided=NA, conf.high.onesided=NA, bandwidth=h,
                     eff.obs=r1$eff.obs,
                     leverage=max(c(0, wt^2/d$w[r1$est_w!=0]^2))/
                         sum(wt^2/d$w[r1$est_w!=0]),
                     cv=NA, alpha=alpha, method=method, M=M[1], M.rf=M[2],
                     M.fs=M[3], first.stage=r1$fs, kernel=kernel_type(kern),
                     p.value=NA)
    structure(list(coefficients=fill_coefs(co), data=d,
                   lm=r1$lm), class="RDResults")
}

## Same for RDTOpt and RD
fill_coefs <- function(co) {
    B <- co$maximum.bias/co$std.error
    cv <- CVb(B, co$alpha)
    co[, c("conf.low", "conf.high", "conf.low.onesided",
           "conf.high.onesided", "cv", "p.value")] <-
        c(co$estimate-cv*co$std.error, co$estimate+cv*co$std.error,
          co$estimate - (B + stats::qnorm(1-co$alpha))*co$std.error,
          co$estimate + (B + stats::qnorm(1-co$alpha))*co$std.error,
          cv, stats::pnorm(B-abs(co$estimate/co$std.error))+
              stats::pnorm(-B-abs(co$estimate/co$std.error)))
    co
}


## Optimal bandwidth selection in nonparametric regression
OptBW <- function(d, M, kern="triangular", opt.criterion, alpha=0.05, beta=0.8,
                  sclass="H", T0=0) {
    ## First check if sigma2 is supplied
    if (is.null(d$sigma2))
        d <- PrelimVar(d, se.initial="EHW")

    ## Objective function for optimizing bandwidth
    obj <- function(h) {
        r <- NPRHonest(d, M, kern, h, alpha=alpha, se.method="supplied.var",
                       sclass=sclass, T0=T0, T0bias=TRUE)$coefficients
        switch(opt.criterion,
               OCI=2*r$maximum.bias+
                   r$std.error * (stats::qnorm(1-alpha)+stats::qnorm(beta)),
               MSE=r$maximum.bias^2+r$std.error^2,
               FLCI=r$conf.high-r$conf.low)
    }
    if (inherits(d, "IP")) {
        hmin <- sort(unique(abs(d$X)))[2]
    } else {
        hmin <- max(unique(d$X[d$p])[2], sort(unique(abs(d$X[d$m])))[2])
    }

    ## Optimize piecewise constant function using modification of golden
    ## search. In fact, the criterion may not be unimodal, so proceed with
    ## caution. (For triangular kernel, it appears unimodal)
    if (kernel_type(kern)=="uniform") {
        supp <- sort(unique(abs(d$X)))
        h <- gss(obj, supp[supp>=hmin])
    } else {
        h <- abs(stats::optimize(obj, interval=c(hmin, max(abs(d$X))),
                                 tol=.Machine$double.eps^0.75)$minimum)
    }

    h
}



#' @export
print.RDResults <- function(x, digits = getOption("digits"), ...) {
    if (!is.null(x$call)) {
        cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
            "\n\n", sep = "")
    }
    fmt <- function(x) format(x, digits=digits, width=digits+1)
    y <- x$coefficients
    cat("Inference for ", y$term, " (using ", y$method,
        " class), confidence level ", 100*(1-y$alpha), "%:\n", sep="")
    nm <- c("Estimate", "Std. Error", "Maximum Bias")
    colnames(y)[2:4] <- nm
    y$"Confidence Interval" <- paste0("(", fmt(y$conf.low), ", ",
                                      fmt(y$conf.high), ")")
    y$OCI <- paste0("(-Inf, ", fmt(y$conf.high.onesided), "), (",
                    fmt(y$conf.low.onesided), ", Inf)")

    print.data.frame(y[, c(nm, "Confidence Interval"), ],
                     digits=digits)
    cat("\nOnesided CIs: ", y$OCI)

    cat("\nNumber of effective observations:", fmt(y$eff.obs))
    par <- paste0(tolower(substr(y$term, 1, 1)), substring(y$term, 2))
    cat("\nMaximal leverage for ", par, ": ", fmt(y$leverage),
        sep="")
    if (!is.null(y$first.stage) && !is.na(y$first.stage))
        cat("\nFirst stage estimate:", fmt(y$first.stage),
            "\nFirst stage smoothness constant M:", fmt(y$M.fs),
            "\nReduced form smoothness constant M:", fmt(y$M.rf))
    else if (!is.null(y[["M"]])) {
        cat("\nSmoothness constant M:", fmt(y$M))
    }

    cat("\nP-value:", fmt(y$p.value), "\n")

    if (!is.null(y$bandwidth)) {
        cat("\nBased on local regression with bandwidth: ", fmt(y$bandwidth),
            ", kernel: ", y$kernel, "\nRegression coefficients:\n", sep="")
        print.default(format(coef(x$lm), digits = max(3L, digits - 3L)),
                      print.gap = 2L, quote = FALSE)
    } else {
        cat("\nSmoothing parameters below and above cutoff: ",
            fmt(y$bandwidth.m), ", ", fmt(y$bandwidth.p), sep="")
    }
    if (inherits(x$na.action, "omit"))
        cat(length(x$na.action), "observations with missing values dropped\n")

    invisible(x)
}
