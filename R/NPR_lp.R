#' Honest inference in RD
#'
#' Basic computing engine called by \code{\link{RDHonest}} to compute honest
#' confidence intervals for local polynomial estimators.
#' @param d object of class \code{"RDData"}
#' @template RDse
#' @template RDoptBW
#' @template RDBW
#' @template RDclass
#' @template Kern
#' @template bwequal
#' @template RDseInitial
#' @return Returns an object of class \code{"RDResults"}, see description in
#'     \code{\link{RDHonest}}
#' @export
NPRHonest.fit <- function(d, M, kern="triangular", h, opt.criterion,
                         bw.equal=TRUE, alpha=0.05, beta=0.8, se.method="nn",
                         J=3, sclass="H", order=1, se.initial="EHW") {
    ## Initial se estimate
    if (is.null(d$sigma2) & (is.null(d$sigma2p) | is.null(d$sigma2m)) &
        ("supplied.var" %in% se.method | missing(h)))
        d <- NPRPrelimVar.fit(d, se.initial=se.initial)

    if (missing(h)) {
        h <- NPROptBW.fit(d, M, kern, opt.criterion, bw.equal, alpha,
                         beta, sclass, order)$h
    } else if (length(h)==1) {
        h <- c(p=unname(h), m=unname(h))
    }

    ## Suppress warnings about too few observations
    r1 <- NPRreg.fit(d, h, kern, order, se.method, TRUE, J)
    if (inherits(d, "LPPData")) {
        w <- r1$w(d$X)
        wt <- w[w!=0]
        xx <- d$X[w!=0]
        nobs <- length(wt)
        ## test for boundary
        bd <- length(unique(d$X>=0)==1)
    } else {
        wp <- r1$wp(d$Xp)
        wm <- r1$wm(d$Xm)
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
        M0 <- if (inherits(d, "FRDData")) (M[1]+M[2]*abs(r1$estimate)) else M
        if(order==0) {
            bias <- Inf
        } else if (sclass=="T")  {
            bias <- M0/2 * (sum(abs(wt*xx^2)))
        } else if (sclass=="H" & order==1 & bd) {
            ## At boundary we know form of least favorable function
            bias <- -M0/2 * (sum(wt*xx^2))
        } else {
            ## Else need to find numerically

            w2p <- function(s) abs(sum((wt*(xx-s))[xx>=s]))
            w2m <- function(s) abs(sum((wt*(s-xx))[xx<=s]))
            bp <- stats::integrate(function(s)
                vapply(s, w2p, numeric(1)), 0, h["p"])$value
            bm <- stats::integrate(function(s)
                vapply(s, w2m, numeric(1)), -h["m"], 0)$value
            bias <- M0*(bp+bm)
        }
        lower <- r1$estimate - bias - stats::qnorm(1-alpha)*sd
        upper <- r1$estimate + bias + stats::qnorm(1-alpha)*sd
        hl <- CVb(bias/sd, alpha)$cv*sd
    }

    ## Finally, calculate coverage of naive CIs
    z <- stats::qnorm(1-alpha/2)
    naive <- stats::pnorm(z-bias/sd)-stats::pnorm(-z- bias/sd)

    structure(list(estimate=r1$estimate, lff=NA, maxbias=bias, sd=sd,
                   lower=lower, upper=upper, hl=hl, eff.obs=r1$eff.obs,
                   hp=unname(h["p"]), hm=unname(h["m"]), naive=naive),
              class="NPRResults")
}


## RDOptBW.fit <- function(d, M, kern="triangular", opt.criterion, bw.equal=TRUE,
##                          alpha=0.05, beta=0.8, sclass="H", order=1,
##                          se.initial="EHW", T0=0) {

##     ## First check if sigma2 is supplied
##     if (is.null(d$sigma2) & (is.null(d$sigma2p) | is.null(d$sigma2m)))
##         d <- NPRPrelimVar.fit(d, se.initial=se.initial)

##     ## Objective function for optimizing bandwidth
##     obj <- function(hp, hm) {
##         r <- NPRHonest.fit(d, M, kern, c(p=abs(hp), m=abs(hm)),
##                           alpha=alpha, se.method="supplied.var",
##                           sclass=sclass, order=order)
##         if (opt.criterion=="OCI") {
##             2*r$maxbias+r$sd*(stats::qnorm(1-alpha)+stats::qnorm(beta))
##         } else if (opt.criterion=="MSE") {
##             r$maxbias^2+r$sd^2
##         } else if (opt.criterion=="FLCI") {
##             r$hl
##         } else {
##             stop(sprintf("optimality criterion %s not yet implemented",
##                          opt.criterion))
##         }
##     }

##     if (bw.equal==FALSE) {
##         r <- stats::optim(c(max(d$Xp)/2, max(abs(d$Xm)/2)),
##                    function(h) obj(h[1], h[2]))
##         if (r$convergence!=0) warning(r$message)
##         hp <- abs(r$par[1])
##         hm <- abs(r$par[2])
##     } else {
##         obj1 <- function(h) obj(h, h)
##         hmin <- max(unique(d$Xp)[order+1], sort(unique(abs(d$Xm)))[order+1])
##         hmax <- max(abs(c(d$Xp, d$Xm)))
##         ## Optimize piecewise constant function using modification of golden
##         ## search. In fact, the criterion may not be unimodal, so proceed with
##         ## caution. (For triangular kernel, it appears unimodal)
##         if (kern=="uniform") {
##             supp <- sort(unique(c(d$Xp, abs(d$Xm))))
##             hp <- hm <- gss(obj1, supp[supp>=hmin])
##         } else {
##             hm <- hp <- abs(stats::optimize(obj1, interval=c(hmin, hmax),
##                                         tol=.Machine$double.eps^0.75)$minimum)
##         }

##     }

##     structure(list(h=c(p=hp, m=hm), sigma2p=d$sigma2p, sigma2m=d$sigma2m),
##               class="NPRBW")
## }


#' Optimal bandwidth selection in RD
#'
#' Basic computing engine called by \code{\link{RDOptBW}} used to find
#' optimal bandwidth
#' @param d object of class \code{"RDData"}
#' @template RDoptBW
#' @template RDclass
#' @template Kern
#' @template bwequal
#' @template RDseInitial
#' @param T0 Inital estimate of beta, for FRD only
#' @return a list with the following elements
#'     \describe{
#'     \item{\code{hp}}{bandwidth for observations above cutoff}
#'
#'     \item{\code{hm}}{bandwidth for observations below cutoff, equal to \code{hp}
#'     unless \code{bw.equal==FALSE}}
#'
#'     \item{\code{sigma2m}, \code{sigma2p}}{estimate of conditional variance
#'      above and below cutoff, from \code{d}}
#'    }
#' @references{
#' \cite{Imbens, Guido, and Kalyanaraman, Karthik,
#' "Optimal bandwidth choice for the regression discontinuity estimator." The
#' Review of Economic Studies 79 (3): 933-959.}
#' }
#' @examples
#' ## Lee data
#' d <- RDData(lee08, cutoff=0)
#' NPROptBW.fit(d, M=0.1, opt.criterion="MSE")[c("hp", "hm")]
#' @export
NPROptBW.fit <- function(d, M, kern="triangular", opt.criterion, bw.equal=TRUE,
                         alpha=0.05, beta=0.8, sclass="H", order=1,
                         se.initial="EHW", T0=0) {

    ## First check if sigma2 is supplied
    if (is.null(d$sigma2) & (is.null(d$sigma2p) | is.null(d$sigma2m)))
        d <- NPRPrelimVar.fit(d, se.initial=se.initial)

    ## Objective function for optimizing bandwidth
    obj <- function(hp, hm) {
        r <- NPRHonest.fit(d, M, kern, c(p=abs(hp), m=abs(hm)),
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

    if (bw.equal==FALSE & !inherits(d, "LPPData")) {
        r <- stats::optim(c(max(d$Xp)/2, max(abs(d$Xm)/2)),
                   function(h) obj(h[1], h[2]))
        if (r$convergence!=0) warning(r$message)
        hp <- abs(r$par[1])
        hm <- abs(r$par[2])
    } else {
        obj1 <- function(h) obj(h, h)
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
            hp <- hm <- gss(obj1, supp[supp>=hmin])
        } else {
            hm <- hp <- abs(stats::optimize(obj1, interval=c(hmin, hmax),
                                        tol=.Machine$double.eps^0.75)$minimum)
        }
    }

    structure(list(h=c(p=hp, m=hm), sigma2p=d$sigma2p, sigma2m=d$sigma2m),
              class="NPRBW")
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

    if (x$hm==x$hp) {
        cat("\n", bw, ": ", format(x$hp, digits=digits), sep="")
    } else {
        cat("\n", bw, " below cutoff: ", format(x$hm, digits=digits), sep="")
        cat("\n", bw, " above cutoff: ", format(x$hp, digits=digits), sep="")
    }
    cat("\nNumber of effective observations:",
        format(x$eff.obs, digits=digits), "\n")

    invisible(x)
}

#' @export
print.NPRBW <- function(x, digits = getOption("digits"), ...) {
    if (!is.null(x$call))
        cat("Call:\n", deparse(x$call), "\n", sep = "", fill=TRUE)
    if (x$h["m"]==x$h["p"]) {
        cat("\nBandwidth: ", format(x$h["m"], digits=digits), "\n", sep="")
    } else {
        cat("\nBandwidth below cutoff: ", format(x$h["m"], digits=digits),
            "\n", sep="")
        cat("\nBandwidth above cutoff: ", format(x$h["p"], digits=digits),
            "\n", sep="")
    }
    invisible(x)
}
