## @param d object of class \code{"RDData"}, \code{"FRDData"}, or
##     \code{"LPPData"}
## @param T0 Initial estimate of the treatment effect for calculating the
##     optimal bandwidth. Only relevant for Fuzzy RD.
## @param T0bias When evaluating the maximum bias of the estimate, use the
##     estimate itself (if \code{T0bias==FALSE}), or use the preliminary
##     estimate \code{T0} (if \code{T0bias==TRUE}). Only relevant for Fuzzy RD.
## @return Returns an object of class \code{"NPRResults"}, see descriptions in
##     \code{\link{RDHonest}}, \code{\link{LPPHonest}}, and
##     \code{\link{FRDHonest}}.
NPRHonest.fit <- function(d, M, kern="triangular", h, opt.criterion, alpha=0.05,
                          beta=0.8, se.method="nn", J=3, sclass="H", T0=0,
                          T0bias=FALSE) {
    if (missing(h))
        h <- NPROptBW.fit(d, M, kern, opt.criterion, alpha, beta, sclass, T0)
    r1 <- NPRreg.fit(d, h, kern, order=1, se.method, J)

    wt <- r1$w[r1$w!=0]
    xx <- d$X[r1$w!=0]
    if (d$class=="IP") {
        nobs <- length(wt)
        ## Are we at a boundary?
        bd <- length(unique(d$X>=0))==1
    } else {
        nobs <- min(sum(xx>=0), sum(xx<0))
        bd <- TRUE
    }

    if (T0bias && d$class=="FRD") {
        ## multiply bias and sd by r1$fs to make if free of first stage
        r1$se <- r1$se*abs(r1$fs)
        M <- unname(c((M[1]+M[2]*abs(T0)), M))
    } else if (!T0bias && d$class=="FRD") {
        M <- unname(c((M[1]+M[2]*abs(r1$estimate)) / abs(r1$fs), M))
    }

    ## Determine bias
    if (nobs==0) {
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
        w2p <- function(s) abs(sum((wt*(xx-s))[xx>=s]))
        w2m <- function(s) abs(sum((wt*(s-xx))[xx<=s]))
        bp <- stats::integrate(function(s) vapply(s, w2p, numeric(1)), 0, h)
        bm <- stats::integrate(function(s) vapply(s, w2m, numeric(1)), -h, 0)
        bias <- M[1]*(bp$value+bm$value)
    }
    lower <- r1$estimate - bias - stats::qnorm(1-alpha)*r1$se
    upper <- r1$estimate + bias + stats::qnorm(1-alpha)*r1$se
    cv <- CVb(bias/r1$se, alpha)

    term <- switch(d$class, IP="Value of conditional mean",
                   SRD="Sharp RD Parameter", "Fuzzy RD Parameter")
    method <- switch(sclass, H="Holder", "Tayor")
    if (d$class!="FRD") M[2:3] <- c(NA, NA)
    d$est_w <- r1$w
    d$sigma2 <- r1$sigma2
    kernel <- if (!is.function(kern)) kern else "user-supplied"
    coef <- data.frame(term=term, estimate=r1$estimate, std.error=r1$se,
                       maximum.bias=bias, conf.low=r1$estimate-cv*r1$se,
                       conf.high=r1$estimate+cv*r1$se, conf.low.onesided=lower,
                       conf.high.onesided=upper, bandwidth=h,
                       eff.obs=r1$eff.obs, leverage=max(r1$w^2)/sum(r1$w^2),
                       cv=cv, alpha=alpha, method=method, M=M[1], M.rf=M[2],
                       M.fs=M[3], first.stage=r1$fs, kernel=kernel)
    structure(list(coefficients=coef, data=d), class="RDResults")
}


## Optimal bandwidth selection in nonparametric regression
NPROptBW.fit <- function(d, M, kern="triangular", opt.criterion, alpha=0.05,
                         beta=0.8, sclass="H", T0=0) {

    ## First check if sigma2 is supplied
    if (is.null(d$sigma2))
        d <- NPRPrelimVar.fit(d, se.initial="EHW")

    ## Objective function for optimizing bandwidth
    obj <- function(h) {
        r <- NPRHonest.fit(d, M, kern, h, alpha=alpha, se.method="supplied.var",
                           sclass=sclass, T0=T0, T0bias=TRUE)$coefficients
        switch(opt.criterion,
               OCI=2*r$maximum.bias+
                   r$std.error*(stats::qnorm(1-alpha)+stats::qnorm(beta)),
               MSE=r$maximum.bias^2+r$std.error^2,
               FLCI=r$conf.high-r$conf.low)
    }
    hmin <- if (d$class=="IP") {
                sort(unique(abs(d$X)))[2]
            } else {
                max(unique(d$X[d$p])[2], sort(unique(abs(d$X[d$m])))[2])
            }

    ## Optimize piecewise constant function using modification of golden
    ## search. In fact, the criterion may not be unimodal, so proceed with
    ## caution. (For triangular kernel, it appears unimodal)
    if (kern=="uniform") {
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
    if (!is.null(x$call))
        cat("Call:\n", deparse(x$call), "\n\n", sep = "", fill=TRUE)
    fmt <- function(x) format(x, digits=digits, width=digits+1)
    y <- x$coefficients
    cat("Estimates (using ", y$method, " class):\n", sep="")
    nm <- c("Parameter", "Estimate", "Std. Error", "Maximum Bias")
    names(y)[1:4] <- nm
    y$"Confidence Interval" <- paste0("(", fmt(y$conf.low), ", ",
                   fmt(y$conf.high), ")")
    y$OCI <- paste0("(-Inf, ", fmt(y$conf.high.onesided), "), (",
                    fmt(y$conf.low.onesided), ", Inf)")
    print.data.frame(y[, c(nm, "Confidence Interval"), ],
                     digits=digits, row.names=FALSE)
    cat("\nOnesided CIs: ", y$OCI)
    if(!is.null(y$bandwidth))
        cat("\nBandwidth: ", fmt(y$bandwidth), ", Kernel: ", y$kernel, sep="")
    else
        cat("\nSmoothing parameters below and above cutoff: ",
            fmt(y$bandwidth.m), ", ",
            fmt(y$bandwidth.m), sep="")

    cat("\nNumber of effective observations:", fmt(y$eff.obs))
    par <- paste0(tolower(substr(y$Parameter, 1, 1)), substring(y$Parameter, 2))
    cat("\nMaximal leverage for ", par, ": ", fmt(y$leverage),
        sep="")
    if(!is.null(y$first.stage) && !is.na(y$first.stage))
        cat("\nFirst stage estimate:", fmt(y$first.stage),
            "\nFirst stage smoothness constant M:", fmt(y$M.fs),
            "\nReduced form smoothness constant M:", fmt(y$M.rf),
            "\n")
    else
        cat("\nSmoothness constant M:", fmt(y$M),
            "\n")

    if (inherits(x$na.action, "omit"))
        cat(length(x$na.action), "observations with missing values dropped\n")



    invisible(x)
}
