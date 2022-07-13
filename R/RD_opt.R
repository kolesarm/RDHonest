## Calculate minimax optimal estimators, tests, and confidence sets for RD
## problem under second-order Taylor smoothness class

## Compute inverse modulus squared divided by 4,
## \eqn{\omega^{-1}(2b)^2/4=sum_{i} g(x_i)^2/sigma^2(x_i)},
## @param d Object of class \code{"RDData"}
## @param f Least favorable function of class \code{"RDLFFunction"}
Q <- function(d, f) sum(f$m(d$Xm)^2/d$sigma2m)+ sum(f$p(d$Xp)^2/d$sigma2p)


## Solution to inverse modulus problem in RD under Taylor(2) class
## Compute function \eqn{g_{b, C}(x)} that solves the inverse modulus problem
## \eqn{omega^{-1}(2b,C)} in RD under second-order Taylor class with smoothness
## parameter \eqn{C} by solving for parameters \eqn{d_{+}}, \eqn{d_{-}}, and
## \eqn{b_{-}}
## @param d Object of class \code{"RDData"}.
## @param b Jump at zero
## @param C smoothness parameter for Taylor class, second derivative at zero is
## bounded by \code{2*C}.
## @return Object of class \code{"RDLFFunction"}
RDgbC <- function(d, b, C) {
    ## Sacks-Ylvisaker function: g(x)=(b+dx-Cx^2)_{+} - (b+dx+Cx^2)_{-}
    SY <- function(x, b, d, C) pmax(b+d*x-C*x^2, 0) + pmin(b+d*x+C*x^2, 0)

    ## Find d for SY function
    dstar <- function(X, b, C, sigma2) {
        FindZero(function(d) sum(SY(X, b, d, C)*X/sigma2))
    }
    ## Find b_{-}
    eq <- function(bm) {
        dm <- dstar(d$Xm, bm,   C, d$sigma2m)
        dp <- dstar(d$Xp, b-bm, C, d$sigma2p)

        sum(SY(d$Xp, b-bm, dp, C) / d$sigma2p) -
            sum(SY(d$Xm, bm, dm, C)/d$sigma2m)
    }

    bm <- FindZero(eq, b)
    bp <- b-bm
    dp <- dstar(d$Xp, bp, C, d$sigma2p)
    dm <- dstar(d$Xm, bm, C, d$sigma2m)

    list(m=function(x) SY(x, bm, dm, C)*(x<=0),
         p=function(x) SY(x, bp, dp, C)*(x>=0))
}


## Solve modulus problem in RD under second-order Taylor class
## Compute function \eqn{g_{b(\delta), C}(x)} that solves the modulus problem
## \eqn{omega(\delta, C)} in RD under second-order Taylor class with smoothness
## parameter \eqn{C}.
## @param delta \eqn{\delta}
## @inheritParams RDgbC
RDLFFunction <- function(d, C, delta) {
    RDgbC(d, FindZero(function(b) 4*Q(d, RDgbC(d, b, C)) - delta^2,
                      negative=FALSE), C)
}

## Compute optimal estimator based on solution to modulus problem, and CIs
## around it
## \eqn{hat{L}_{delta, C}}
## @param d Object of class \code{"RDData"}
## @param f RDLFFunction
RDTEstimator <- function(d, f, alpha, se.method, J) {
    den <- sum(f$p(d$Xp) / d$sigma2p) # denominator
    ## By (S3) in supplement of 1511.06028v2, this = sum(f$m(d$Xm) / d$sigma2m)

    Wp <- f$p(d$Xp) / (d$sigma2p*den)
    Wm <- f$m(d$Xm) / (d$sigma2m*den)

    q <- Q(d, f)
    ## If not NN then it's suppled var
    if (se.method=="nn") {
        d$sigma2p <- sigmaNN(d$Xp, d$Yp, J=J)
        d$sigma2m <- sigmaNN(d$Xm, d$Ym, J=J)
    }
    ## Drop observations with zero weight
    d$Xp <- d$Xp[Wp!=0]
    d$Xm <- d$Xm[Wm!=0]
    d$Yp <- d$Yp[Wp!=0]
    d$Ym <- d$Ym[Wm!=0]
    d$sigma2p <- d$sigma2p[Wp!=0]
    d$sigma2m <- d$sigma2m[Wm!=0]
    Wp <- Wp[Wp!=0]
    Wm <- Wm[Wm!=0]
    Lhat <- sum(Wp * d$Yp) - sum(Wm * d$Ym)

    b <- f$m(0)+f$p(0)

    sd <-  sqrt(sum(Wp^2 * d$sigma2p) + sum(Wm^2 * d$sigma2m))
    maxbias <- b - q/den                # b-q/den
    lower <- Lhat - maxbias - stats::qnorm(1-alpha)*sd
    upper <- Lhat + maxbias + stats::qnorm(1-alpha)*sd
    hl <- CVb(maxbias/sd, alpha) * sd # Half-length

    ## Effective number of observations
    ## eff.obs <- 1/sum(Wp^2) + 1/sum(Wm^2)
    coef <- data.frame(term="Sharp RD parameter", estimate=Lhat, std.error=sd,
                       maximum.bias=maxbias, conf.low=Lhat-hl,
                       conf.high=Lhat+hl, conf.low.onesided=lower,
                       conf.high.onesided=upper, eff.obs=NA, #TODO
                       cv=CVb(maxbias/sd, alpha), alpha=alpha, method="Taylor")

    structure(list(coefficients=coef, data=d, delta=sqrt(4*q), omega=2*b),
              class="RDResults")
}

## Optimal inference in RD under Taylor class
RDTOpt.fit <- function(d, M, opt.criterion, alpha, beta, se.method, J) {
    ## First check if sigma2 is supplied
    if (is.null(d$sigma2p) || is.null(d$sigma2m))
        d <- NPRPrelimVar.fit(d, se.initial="EHW")

    C <-  M/2
    ## Find optimal delta, see Supplement to 1511.06028v2
    if (opt.criterion=="OCI") {
        lff <- RDLFFunction(d, C, stats::qnorm(1-alpha)+stats::qnorm(beta))
    } else if (opt.criterion=="MSE") {
        eq <- function(b) {
            r <- RDgbC(d, b, C)
            C*sum((d$Xp)^2*abs(r$p(d$Xp)) / d$sigma2p) +
                C*sum(d$Xm^2 * abs(r$m(d$Xm)) / d$sigma2m) - 1
        }
        lff <- RDgbC(d, FindZero(eq, negative=FALSE), C)
    } else if (opt.criterion=="FLCI") {
        eq <- function(b) {
            r <- RDgbC(d, b, C)
            q1 <- Q(d, r)
            ## Instead of using the above formula from page S5 of 1511.06028v2,
            ## optimize half-length directly
            den <- sum(r$p(d$Xp) / d$sigma2p) # denominator
            hse <- sqrt(q1)/den                   # standard deviation
            maxbias <- b - q1/den                # b-q/den
            CVb(maxbias/hse, alpha) * hse # Half-length
        }
        ## eq is convex, start around MSE optimal b
        bs <- RDTOpt.fit(d, M, "MSE", alpha, beta, se.method, J)$omega/2
        lff <- RDgbC(d, stats::optimize(eq, c(bs/2, 3*bs/2))$minimum, C)
    }

    ## Compute optimal estimator
    r <- RDTEstimator(d, lff, alpha, se.method, J)
    ## Two bandwidths in this case
    r$coefficients$bandwidth.m <- (lff$m(0)/C)^(1/2)
    r$coefficients$bandwidth.p <- (lff$p(0)/C)^(1/2)
    r
}


#' Finite-sample efficiency bounds for minimax CIs
#'
#' Compute efficiency of minimax one-sided CIs at constant functions, or
#' efficiency of two-sided fixed-length CIs at constant functions under
#' second-order Taylor smoothness class.
#'
#' @param object An object of class \code{"RDResults"}, typically a result of a
#'     call to \code{\link{RDHonest}}.
#' @param opt.criterion Either \code{"FLCI"} for computing efficiency of
#'     two-sided CIs, or else \code{"OCI"} for minimax one-sided CIs.
#' @param beta Determines quantile of excess length for evaluating minimax
#'     efficiency of one-sided CIs. Ignored if \code{opt.criterion=="FLCI"}.
#' @return Efficiency bound, a numeric vector of length one.
#' @references{
#'
#' \cite{Timothy B. Armstrong and Michal Kolesár. Optimal inference in a class
#' of regression models. Econometrica, 86(2):655–683, March 2018.
#' \doi{10.3982/ECTA14434}}
#'
#' }
#' @examples
#' r <- RDHonest(voteshare ~ margin, data=lee08, M=0.1, h=2)
#' RDTEfficiencyBound(r, opt.criterion="OCI")
#' @export
RDTEfficiencyBound <- function(object, opt.criterion="FLCI", beta=0.5) {
    d <- object$data
    alpha <- object$coefficients$alpha
    C <- object$coefficients$M/2
    d <- NPRPrelimVar.fit(d, se.initial="EHW")

    if (opt.criterion=="OCI") {
        delta <- stats::qnorm(1-alpha)+stats::qnorm(beta)
        r1 <- RDTEstimator(d, RDLFFunction(d, C, delta), alpha,
                           se.method="supplied.var")
        r2 <- RDTEstimator(d, RDLFFunction(d, C, 2*delta), alpha,
                           se.method="supplied.var")
        return(r2$omega/(r1$delta*r1$coefficients$std.error+r1$omega))
    } else {
        ## From proof of Pratt result, it follows that the expected length is
        ## int pnorm(z_{1-alpha}-delta_t) dt, where delta_t is value of inverse
        ## two-class modulus omega{delta, F, 1}, or omega{delta, 1, F},
        ## depending on whether t > 0 or < 0.
        deltat <- function(t) sqrt(Q(d, RDgbC(d, t, C))) # delta_t
        integrand <- function(t) {
            stats::pnorm(stats::qnorm(1-alpha)-vapply(t, deltat, numeric(1)))
        }
        ## By symmetry, half-length is given by value of integral over R_+. The
        ## integrand equals 1-alpha at zero, need upper cutoff
        upper <- 10
        while(integrand(upper)>1e-10) upper <- 2*upper
        den <- RDTOpt.fit(d, 2*C, opt.criterion="FLCI", alpha, beta,
                          se.method="supplied.var")$coefficients
        den <- (den$conf.high-den$conf.low)/2
        return(stats::integrate(integrand, 1e-6, upper)$value / den)
    }
}
