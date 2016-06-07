#' Honest inference in RD based on local polynomial estimator
#'
#' Calculate estimators and one- and two-sided CIs based on local polynomial
#' estimator in RD under second-order Taylor or Hölder smoothness class. If
#' bandwidth is not supplied it is calculated to be optimal for a given
#' performance criterion, as specified by \code{opt.criterion} using the
#' function \code{\link{RDOptBW}}.
#'
#' @template RDFormula
#' @template RDse
#' @template RDoptBW
#' @template RDBW
#' @template RDclass
#' @template Kern
#' @return Returns an object of class \code{"RDResults"}. The function \code{print}
#'     can be used to obtain and print a summary of the results. An object of
#'     class \code{"RDResults"} is a list containing the following components (TODO)
#'
#'     \describe{
#'   \item{\code{estimate}}{Point estimate. This estimate is MSE-optimal if
#'                   \code{what="Estimation"}}
#'
#'   \item{\code{lower}, \code{upper}}{Lower (upper) end-point of a one-sided CI
#'         based on \code{estimate}. This CI is optimal if \code{what=="OCI"}}
#'
#'   \item{\code{hl}}{Half-length of a two-sided CI based on \code{estimate},
#'             so that the
#'             CI is given by \code{c(estimate-hl, estimate+hl)}. The CI is optimal
#'             if \code{what="FLCI"}}
#'   \item{\code{maxbias}}{Maximum bias of \code{estimate}}
#'   \item{\code{sd}}{Standard deviation of estimate}
#'   \item{\code{eff.obs}}{Effective number of observations used by \code{estimate}}
#'   \item{\code{hp}, \code{hm}}{Bandwidths used}
#'   \item{\code{naive}}{Coverage of naive CI that uses \code{qnorm(1-alpha/2)} as
#'                critical value}
#' }
#' @seealso \code{\link{RDOptBW}}
#' @examples
#'
#' ## TODO
#' @export
RDHonest <- function(formula, data, subset, cutoff=0, M, kern="triangular",
                     sigma2, na.action, opt.criterion, bw.equal=TRUE, hp, hm=hp,
                     se.method="nn", alpha=0.05, beta=0.8, J=3, sclass="H",
                     order=1) {

    ## construct model frame
    cl <- mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "sigma2"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    d <- RDData(mf, cutoff)

    if (!missing(hp)) {
        ret <- RDHonest.fit(d, M, kern, hp, hm, alpha=alpha,
                            se.method=se.method, J=J, sclass=sclass,
                            order=order)
    } else {
        ret <- RDHonest.fit(d, M, kern, opt.criterion=opt.criterion,
                            bw.equal=bw.equal, alpha=alpha, beta=beta,
                            se.method=se.method, J=J,
                            sclass=sclass, order=order)
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
#' @examples
#'
#' ## TODO
#' @export
RDOptBW <- function(formula, data, subset, cutoff=0, M, kern="triangular",
                    sigma2, na.action, opt.criterion, bw.equal=TRUE,
                    alpha=0.05, beta=0.8, sclass="H", order=1) {

    ## construct model frame
    cl <- mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "sigma2"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    d <- RDData(mf, cutoff)

    ret <- RDOptBW.fit(d, M, kern, opt.criterion, bw.equal, alpha, beta, sclass,
                       order)
    class(ret) <- "RDBW"
    ret$call <- cl
    ret$na.action <- attr(mf, "na.action")

    ret
}


#' Basic computing engine called by \code{\link{RDHonest}} to compute honest
#' confidence intervals
#' @param d object of class \code{"RDData"}
#' @template RDse
#' @template RDoptBW
#' @template RDBW
#' @template RDclass
#' @template Kern
#' @return TODO
#' @examples
#'
#' ## TODO
#' @importFrom stats pnorm qnorm
#' @export
RDHonest.fit <- function(d, M, kern="triangular", hp, hm=hp, opt.criterion,
                         bw.equal=TRUE, alpha=0.05, beta=0.8, se.method="nn",
                         J=3, sclass="H", order=1) {
    CheckClass(d, "RDData")

    if (missing(hp)) {
        r <- RDOptBW.fit(d, M, kern, opt.criterion, bw.equal, alpha,
                         beta, sclass, order)
        hp <- r$hp
        hm <- r$hm
    }

    ## Suppress warnings about too few observations
    r1 <- RDLPreg(d, hp, kern, order, hm, se.method, TRUE, J)
    wp <- r1$wp(d$Xp)
    wm <- r1$wm(d$Xm)

    ## If bandwidths too small
    if (sum(wp>0)==0 | sum(wm>0)==0) {
        ## big bias / sd
        bias <- sd <- upper <- hl <- sqrt(.Machine$double.xmax/10)
        lower <- -upper
    } else {
        sd <- r1$se[se.method]
        if(sclass=="T")  {
            bias <- M/2 * (sum(abs(wp*d$Xp^2)) + sum(abs(wm*d$Xm^2)))
        } else if (sclass=="H" & order==1) {
            bias <- -M/2 * (sum(wp*d$Xp^2) + sum(wm*d$Xm^2))
        } else if (sclass=="H" & order==2) {
            ## need to find numerically
            bp <- function(b) -2*sum(wp*(pmax(d$Xp-b, 0)^2))
            bm <- function(b) -2*sum(wm*(pmax(abs(d$Xm)-b, 0)^2))
            bias <- -(M/2)*(CarefulOptim(bp, c(0, hp), k=5)$objective +
                            CarefulOptim(bm, c(0, hm), k=5)$objective)
        } else {
            stop("Don't know how to compute bias for
                  specified sclass and order.")
        }

        lower <- r1$estimate - bias - qnorm(1-alpha)*sd
        upper <- r1$estimate + bias + qnorm(1-alpha)*sd
        hl <- CVb(bias/sd, alpha)$cv*sd
    }

    ## Finally, calculate coverage of naive CIs
    z <- qnorm(1-alpha/2)
    naive <- pnorm(z-bias/sd)-pnorm(-z- bias/sd)

    structure(list(estimate=r1$estimate, lff=NA, maxbias=bias, sd=sd,
                   lower=lower, upper=upper, hl=hl, eff.obs=r1$eff.obs,
                   hp=hp, hm=hm, naive=naive),
              class="RDResults")
}


#' Basic computing engine called by \code{\link{RDOptBW}} used to find
#' optimal bandwidth
#' @param d object of class \code{"RDData"}
#' @template RDoptBW
#' @template RDclass
#' @template Kern
#' @return TODO
#' @examples
#' ## TODO
#' @export
RDOptBW.fit <- function(d, M, kern="triangular", opt.criterion,
                        bw.equal=TRUE, alpha=0.05, beta=0.8,
                        sclass="H", order=1) {

    ## First check if sigma2 is supplied, if not, use nn method
    if (is.null(d$sigma2p) | is.null(d$sigma2m)) {
        ## Silverman pilot bandwidth for uniform kernel
        ## TODO: something more reasonable
        X <- c(d$Xm, d$Xp)
        h1 <- 1.84*stats::sd(X)/sum(length(X))^(1/5)
        ## Make sure this results in enough distinct values on either side of
        ## threshold

        h1 <- max(h1,
                  sort(unique(d$Xp))[order+1],
                  abs(sort(unique(d$Xm), decreasing=TRUE)[order+1]))

        r1 <- RDLPreg(d=d, hp=h1, kern=kern, order=1, se.method="nn")
        d$sigma2p <- rep(mean(r1$sigma2p), length(d$Xp))
        d$sigma2m <- rep(mean(r1$sigma2m), length(d$Xm))
    }

    ## Objective function for optimizing bandwidth
    obj <- function(hp, hm) {
        r <- RDHonest.fit(d, M, kern, abs(hp), abs(hm),
                          alpha=alpha, se.method="supplied.var",
                          sclass=sclass, order=order)
        if (opt.criterion=="OCI") {
            2*r$maxbias+r$sd*(qnorm(1-alpha)+qnorm(beta))
        } else if (opt.criterion=="MSE") {
            r$maxbias^2+r$sd^2
        } else if (opt.criterion=="FLCI") {
            r$hl
        } else {
            stop(sprintf("optimality criterion %s not yet implemented",
                         opt.criterion))
        }
    }

    if (bw.equal==FALSE) {
        r <- stats::optim(c(max(d$Xp)/2, max(abs(d$Xm)/2)),
                   function(h) obj(h[1], h[2]))
        if (r$convergence!=0) warning(r$message)
        hp <- abs(r$par[1])
        hm <- abs(r$par[2])
    } else {
        hm <- hp <- abs(stats::optimize(function(h)
            obj(h, h), interval=c(0, max(abs(c(d$Xp, d$Xm)))),
                                        tol=tol)$minimum)
    }

    list(hp=hp, hm=hm, sigma2p=d$sigma2p, sigma2m=d$sigma2m)
}

#' Imbens and Kalyanaraman bandwidth
#'
#' Calculate bandwidth for sharp RD based on local linear regression using
#' method by Imbens and Kalyanaraman (2012, ReStud)
#' @param d object of class \code{"RDData"}
#' @template Kern
#' @param verbose Print details of calculation?
#' @return IK bandwidth
#' @importFrom stats coef lm
#' @export
IKBW.fit <- function(d, kern="triangular", order=1, verbose=FALSE) {
    if (order!=1)
        stop("Only works for local linear regression.")

    X <- c(d$Xm, d$Xp)
    Nm <- length(d$Xm)
    Np <- length(d$Xp)
    N <- Nm+Np

    ## STEP 0: Kernel constant
    if (is.character(kern)) {
        s <- kernC[with(kernC,
                        order==1 & boundary==TRUE & kernel==kern, ), ]
    } else if (is.function(kern)) {
        ke <- EqKern(kern, boundary=TRUE, order=1)
        s <- list(mu2=KernMoment(ke, moment=2, boundary=TRUE, "raw"),
                  nu0=KernMoment(ke, moment=0, boundary=TRUE, "raw2"))
    }
    const <- (s$nu0/s$mu2^2)^(1/5)

    ## STEP 1: Estimate f(0), sigma^2_(0) and sigma^2_+(0)
    h1 <- 1.84*stats::sd(X)/N^(1/5)  ## Silverman pilot bandwidth for uniform kernel
    f0 <- sum(X >= -h1 & X<= h1) / (2*N*h1)
    varm <- stats::var(d$Ym[d$Xm >= -h1]) # Equation (12) in IK
    varp <- stats::var(d$Yp[d$Xp <= h1])

    ## STEP 2: Estimate second derivatives m_{+}^(2) and m_{-}^(2)

    ## Estimate third derivative using 3rd order polynomial: Equation (14)
    m3 <- 6*coef(lm(I(c(d$Ym, d$Yp)) ~ I(X>=0) + X + I(X^2) + I(X^3)))[5]

    ## Left and right bandwidths, Equation (15) and page 956.
    ## Optimal constant based on one-sided uniform Kernel, 7200^(1/7),
    h2m <- 7200^(1/7) * (varm/(f0*m3^2))^(1/7) * Nm^(-1/7)
    h2p <- 7200^(1/7) * (varp/(f0*m3^2))^(1/7) * Np^(-1/7)

    ## estimate second derivatives by local quadratic
    m2m <- 2*coef(lm(d$Ym ~ d$Xm + I(d$Xm^2), subset=(d$Xm >= -h2m)))[3]
    m2p <- 2*coef(lm(d$Yp ~ d$Xp + I(d$Xp^2), subset=(d$Xp <= h2p)))[3]

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




#' @export
print.RDBW <- function(x, digits = getOption("digits"), ...) {
    cat("Call:\n", deparse(x$call), "\n\n", sep = "")
    cat("Bandwidth below cutoff: ", format(x$hm, digits=digits))
    cat("\nBandwidth above cutoff: ", format(x$hp, digits=digits))
    if (x$hm==x$hp) {
        cat(" (Bandwidths are the same)\n\n")
    } else {
        cat(" (Bandwidths are different)\n\n")
    }
    invisible(x)
}


#' @export
print.RDResults <- function(x, digits = getOption("digits"), ...) {
    if (!is.null(x$call))
        cat("Call:\n", deparse(x$call), "\n\n", sep = "")

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

    cat("\nBandwidth below cutoff: ", format(x$hm, digits=digits))
    cat("\nBandwidth above cutoff: ", format(x$hp, digits=digits))
    if (x$hm==x$hp) {
        cat(" (Bandwidths are the same)\n")
    } else {
        cat(" (Bandwidths are different)\n")
    }
    cat("Number of effective observations:",
        format(x$eff.obs, digits=digits), "\n")

    invisible(x)
}
