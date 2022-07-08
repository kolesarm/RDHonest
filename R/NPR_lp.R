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

    ## Suppress warnings about too few observations
    r1 <- NPRreg.fit(d, h, kern, order=1, se.method, TRUE, J)
    if (inherits(d, "LPPData")) {
        w <- r1$w
        wt <- w[w!=0]
        xx <- d$X[w!=0]
        nobs <- length(wt)
        ## Are we at a boundary?
        bd <- length(unique(d$X>=0))==1
    } else {
        wp <- r1$wp
        wm <- r1$wm
        wt <- c(wm[wm!=0], wp[wp!=0])
        xx <-  c(d$Xm[wm!=0], d$Xp[wp!=0])
        nobs <- min(sum(wp!=0), sum(wm!=0))
        bd <- TRUE
    }

    ## If bandwidths too small
    if (nobs==0) {
        ## big bias / sd
        bias <- sd <- upper <- hl <- sqrt(.Machine$double.xmax/10)
        lower <- -upper
    } else {
        sd <- r1$se[se.method]
        if (T0bias==TRUE && inherits(d, "FRDData")) {
            ## rescale bias and sd to make it free of r1$fs, use T0
            sd <- sd*abs(r1$fs)
            M <- unname((M[1]+M[2]*abs(T0)))
        } else if (T0bias==FALSE && inherits(d, "FRDData")) {
            M <- unname(M[1]+M[2]*abs(r1$estimate)) / abs(r1$fs)
        }

        if (sclass=="T")  {
            bias <- M/2 * (sum(abs(wt*xx^2)))
        } else if (sclass=="H" && bd) {
            ## At boundary we know form of least favorable function
            bias <- -M/2 * (sum(wt*xx^2))
        } else {
            ## Else need to find numerically (same formula if order=2)
            w2p <- function(s) abs(sum((wt*(xx-s))[xx>=s]))
            w2m <- function(s) abs(sum((wt*(s-xx))[xx<=s]))
            bp <- stats::integrate(function(s) vapply(s, w2p, numeric(1)), 0,
                                   h)$value
            bm <- stats::integrate(function(s) vapply(s, w2m, numeric(1)), -h,
                                   0)$value
            bias <- M*(bp+bm)
        }
        lower <- r1$estimate - bias - stats::qnorm(1-alpha)*sd
        upper <- r1$estimate + bias + stats::qnorm(1-alpha)*sd
        hl <- CVb(bias/sd, alpha)*sd
    }
    term <- if(inherits(d, "LPPData")) {
                "Value of conditional mean"
            } else if (inherits(d, "RDData")) {
                "Sharp RD Parameter"
            } else {
                "Fuzzy RD Parameter"
            }
    coef <- data.frame(
        term=term,
        estimate=r1$estimate,
        std.error=sd,
        maximum.bias=bias,
        conf.low=r1$estimate-hl,
        conf.high=r1$estimate+hl,
        conf.low.onesided=lower,
        conf.high.onesided=upper,
        bandwidth=h,
        eff.obs=r1$eff.obs, # TODO
        cv=NA,
        alpha=alpha,
        method=if (sclass=="H") "Holder" else "Taylor",
        M=M
    )

    structure(list(coefficients=coef,
                   fs=r1$fs), class="RDResults")
}


## Optimal bandwidth selection in nonparametric regression
NPROptBW.fit <- function(d, M, kern="triangular", opt.criterion, alpha=0.05,
                         beta=0.8, sclass="H", T0=0) {

    ## First check if sigma2 is supplied
    if (is.null(d$sigma2) && (is.null(d$sigma2p) || is.null(d$sigma2m)))
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

    if (inherits(d, "LPPData")) {
        hmin <- sort(unique(abs(d$X)))[2]
        X <- d$X
    } else {
        hmin <- max(unique(d$Xp)[2], sort(unique(abs(d$Xm)))[2])
        X <- c(d$Xp, d$Xm)
    }

    ## Optimize piecewise constant function using modification of golden
    ## search. In fact, the criterion may not be unimodal, so proceed with
    ## caution. (For triangular kernel, it appears unimodal)
    if (kern=="uniform") {
        supp <- sort(unique(abs(X)))
        h <- gss(obj, supp[supp>=hmin])
    } else {
        h <- abs(stats::optimize(obj, interval=c(hmin, max(abs(X))),
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
    cat("\nBandwidth: ", format(y$bandwidth, digits=digits), sep="")
    cat("\nNumber of effective observations:",
        format(y$eff.obs, digits=digits), "\n")
    if (inherits(x$na.action, "omit"))
        cat(length(x$na.action), "observations with missing values dropped\n")

    invisible(x)
}
