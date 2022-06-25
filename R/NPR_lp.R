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
                          beta=0.8, se.method="nn", J=3, sclass="H", order=1,
                          se.initial="EHW", T0=0, T0bias=FALSE) {
    ## Initial se estimate
    if (is.null(d$sigma2) & (is.null(d$sigma2p) | is.null(d$sigma2m)) &
        ("supplied.var" %in% se.method | missing(h)))
        d <- NPRPrelimVar.fit(d, se.initial=se.initial)

    if (missing(h))
        h <- NPROptBW.fit(d, M, kern, opt.criterion, alpha, beta, sclass, order,
                          se.initial, T0)$h

    ## Suppress warnings about too few observations
    r1 <- NPRreg.fit(d, h, kern, order, se.method, TRUE, J)
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
        if (T0bias==TRUE & inherits(d, "FRDData")) {
            ## rescale bias and sd to make it free of r1$fs, use T0
            sd <- sd*abs(r1$fs)
            M <- unname((M[1]+M[2]*abs(T0)))
        } else if (T0bias==FALSE & inherits(d, "FRDData")) {
            M <- unname(M[1]+M[2]*abs(r1$estimate)) / abs(r1$fs)
        }

        if(order==0) {
            bias <- Inf
        } else if (sclass=="T")  {
            bias <- M/2 * (sum(abs(wt*xx^2)))
        } else if (sclass=="H" & order==1 & bd) {
            ## At boundary we know form of least favorable function
            bias <- -M/2 * (sum(wt*xx^2))
        } else {
            ## Else need to find numerically
            w2p <- function(s) abs(sum((wt*(xx-s))[xx>=s]))
            w2m <- function(s) abs(sum((wt*(s-xx))[xx<=s]))
            bp <- stats::integrate(function(s)
                vapply(s, w2p, numeric(1)), 0, h)$value
            bm <- stats::integrate(function(s)
                vapply(s, w2m, numeric(1)), -h, 0)$value
            bias <- M*(bp+bm)
        }
        lower <- r1$estimate - bias - stats::qnorm(1-alpha)*sd
        upper <- r1$estimate + bias + stats::qnorm(1-alpha)*sd
        hl <- CVb(bias/sd, alpha)*sd
    }

    ## Finally, calculate coverage of naive CIs
    z <- stats::qnorm(1-alpha/2)
    naive <- stats::pnorm(z-bias/sd)-stats::pnorm(-z- bias/sd)

    structure(list(estimate=r1$estimate, maxbias=bias, sd=sd, lower=lower,
                   upper=upper, hl=hl, eff.obs=r1$eff.obs, h=h, naive=naive,
                   fs=r1$fs), class="NPRResults")
}


## Optimal bandwidth selection in nonparametric regression
##
## Basic computing engine to compute the optimal bandwidth
NPROptBW.fit <- function(d, M, kern="triangular", opt.criterion, alpha=0.05,
                         beta=0.8, sclass="H", order=1, se.initial="EHW", T0=0) {

    ## First check if sigma2 is supplied
    if (is.null(d$sigma2) & (is.null(d$sigma2p) | is.null(d$sigma2m)))
        d <- NPRPrelimVar.fit(d, se.initial=se.initial)

    ## Objective function for optimizing bandwidth
    obj <- function(h) {
        r <- NPRHonest.fit(d, M, kern, h, alpha=alpha, se.method="supplied.var",
                           sclass=sclass, order=order, T0=T0, T0bias=TRUE)
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

    if (inherits(d, "LPPData")) {
        hmin <- sort(unique(abs(d$X)))[order+1]
        hmax <- max(abs(d$X))
    } else {
        hmin <- max(unique(d$Xp)[order+1], sort(unique(abs(d$Xm)))[order+1])
        hmax <- max(abs(c(d$Xp, d$Xm)))
    }
    ## Optimize piecewise constant function using modification of golden
    ## search. In fact, the criterion may not be unimodal, so proceed with
    ## caution. (For triangular kernel, it appears unimodal)
    if (kern=="uniform") {
        supp <- if (inherits(d, "LPPData")) {
                    sort(unique(abs(d$X)))
                } else {
                    sort(unique(c(d$Xp, abs(d$Xm))))
                }
        h <- gss(obj, supp[supp>=hmin])
    } else {
        h <- abs(stats::optimize(obj, interval=c(hmin, hmax),
                                 tol=.Machine$double.eps^0.75)$minimum)
    }

    list(h=h, sigma2p=d$sigma2p, sigma2m=d$sigma2m, sigma2=d$sigma2)
}



#' @export
print.NPRResults <- function(x, digits = getOption("digits"), ...) {
    if (!is.null(x$call))
        cat("Call:\n", deparse(x$call), "\n\n", sep = "", fill=TRUE)
    bw <- if (class(x$lff) != "RDLFFunction")
              "Bandwidth" else "Smoothing parameter"

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

    cat("\n", bw, ": ", format(x$h, digits=digits), sep="")

    cat("\nNumber of effective observations:",
        format(x$eff.obs, digits=digits), "\n")

    invisible(x)
}
