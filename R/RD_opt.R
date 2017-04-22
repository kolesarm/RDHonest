## Calculate minimax optimal estimators, tests, and confidence sets for RD
## problem under second-order Taylor smoothness class

#' Compute inverse modulus squared divided by 4
#'
#' Computes \eqn{\omega^{-1}(2b)^2/4=sum_{i} g(x_i)^2/sigma^2(x_i)},
#' @param d RDData
#' @param f Least favorable function of class "RDLFFunction"
#' @keywords internal
Q <- function(d, f) {
    if (is.null(d$sigma2p) | is.null(d$sigma2m))
        stop("variance function not supplied")

    sum(f$m(d$Xm)^2/d$sigma2m)+ sum(f$p(d$Xp)^2/d$sigma2p)
}

#' Solution to inverse modulus problem in RD under Taylor(2) class
#'
#' Compute function \eqn{g_{b, C}(x)} that solves the inverse modulus problem
#' \eqn{omega^{-1}(2b,C)} in RD under second-order Taylor class with smoothness
#' parameter \eqn{C} by solving for parameters \eqn{d_{+}}, \eqn{d_{-}}, and
#' \eqn{b_{-}}
#' @param d Object of class RDData.
#' @param b Jump at zero
#' @param C smoothness parameter for Taylor class, second derivative at zero is
#' bounded by \code{2*C}.
#' @return Object of class \code{"RDLFFunction"}
#' @keywords internal
RDgbC <- function(d, b, C) {
    ## Sacks-Ylvisaker function: g(x)=(b+dx-Cx^2)_{+} - (b+dx+Cx^2)_{-}
    SY <- function(x, b, d, C) pmax(b+d*x-C*x^2, 0) + pmin(b+d*x+C*x^2, 0)

    ## Find d for SY function
    dstar <- function(X, b, C, sigma2)
        FindZero(function(d) sum(SY(X, b, d, C)*X/sigma2))

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
    cs <- c(bm, bp, dm, dp)
    names(cs) <- c("b.minus", "b.plus", "d.minus", "d.plus")

    structure(list(m=function(x) SY(x, bm, dm, C)*(x<=0),
                   p=function(x) SY(x, bp, dp, C)*(x>=0), c=cs),
              class="RDLFFunction")
}


#' Solve modulus problem in RD under second-order Taylor class
#'
#' Compute function \eqn{g_{b(\delta), C}(x)} that solves the modulus problem
#' \eqn{omega(\delta, C)} in RD under second-order Taylor class with smoothness
#' parameter \eqn{C}.
#' @param delta \eqn{\delta}
#' @inheritParams RDgbC
#' @return Object of class \code{"RDLFFunction"}
#' @keywords internal
RDLFFunction <- function(d, C, delta)
    RDgbC(d, FindZero(function(b) 4*Q(d, RDgbC(d, b, C)) - delta^2,
                      negative=FALSE), C)


#' Compute optimal estimator based on solution to modulus problem, and CIs
#' around it
#'
#' \eqn{hat{L}_{delta, C}}
#' @param d RDData
#' @param f RDLFFunction
#' @keywords internal
RDEstimator <- function(d, f, alpha=0.05) {
    den <- sum(f$p(d$Xp) / d$sigma2p) # denominator
    Lhat <- (sum(f$p(d$Xp) * d$Yp / d$sigma2p) -
                 sum(f$m(d$Xm) * d$Ym / d$sigma2m)) / den
    q <- Q(d, f)
    b <- f$m(0)+f$p(0)
    sd <- sqrt(q)/den                   # standard deviation
    maxbias <- b - q/den                # b-q/den
    lower <- Lhat - maxbias - qnorm(1-alpha)*sd
    upper <- Lhat + maxbias + qnorm(1-alpha)*sd
    hl <- RDHonest::CVb(maxbias/sd, alpha)$cv * sd # Half-length

    ## Effective number of observations
    ## 1/\sum_i w_i^2, since var(Lhat)=sigma^2 *\sum_i w_i^2 under homo
    eff.obs <- 1/(sum( (f$p(d$Xp)/(d$sigma2p*den))^2) +
        sum( (f$m(d$Xm)/(d$sigma2m*den))^2))

    structure(list(estimate=Lhat, lff=f, maxbias=maxbias, sd=sd, lower=lower,
                   upper=upper, hl=hl, delta=sqrt(4*q), omega=2*b,
                   eff.obs=eff.obs),
              class="RDResults")
}

#' Finite-sample optimal inference in RD under second-order Taylor class
#'
#' Compute optimal CI or optimal estimator (depending on value of
#' \code{opt.criterion}) in RD under second-order Taylor class. Variance of the
#' estimator is computed using the conditional variance supplied by \code{d}
#'
#' @param d object of class \code{"RDData"}
#' @template RDoptBW
#' @param M Bound on second derivative of the conditional mean function.
#' @return Returns an object of class \code{"RDResults"}, see description in
#'     \code{\link{RDHonest}}
#' @export
RDTOpt.fit <- function(d, M, opt.criterion, alpha=0.05, beta=0.5) {
    ## We just need to find optimal delta, see end of appendix C for formulas
    C <-  M/2
    if (opt.criterion=="OCI") {
        lff <- RDLFFunction(d, C, qnorm(1-alpha)+qnorm(beta))
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
            q <- Q(d, r)
            ## b*bnmflci(sqrt(q), alpha=alpha)$c *
            ##     sum(RDgbC(d, b, C)$p(d$Xp) / d$sigma2p) - q

            ## Instead of using the above formula from page S5 of 1511.06028v2,
            ## optimize half-length directly
            den <- sum(r$p(d$Xp) / d$sigma2p) # denominator
            sd <- sqrt(q)/den                   # standard deviation
            maxbias <- b - q/den                # b-q/den
            CVb(maxbias/sd, alpha)$cv * sd # Half-length
        }
        ## find interval to optimize over
        b.start <- RDTOpt.fit(d, M, opt.criterion="MSE")$omega/2
        lff <- RDgbC(d, CarefulOptim(eq, c(b.start/2,2*b.start))$minimum, C)
    }

    r <- RDEstimator(d, lff, alpha=alpha)
    r$hm <- (r$lff$m(0)/C)^(1/2)
    r$hp <- (r$lff$p(0)/C)^(1/2)
    r
}


#' Finite-sample efficiency bounds for minimax CIs
#'
#' Compute efficiency of minimax one-sided CIs at constant functions and
#' half-length of Pratt CIs.
#'
#' @param d object of class \code{"RDData"}
#' @template RDoptBW
#' @param M Bound on second derivative of the conditional mean function.
#' @export
RDTEfficiencyBound <- function(d, M, opt.criterion="FLCI",
                               alpha=0.05, beta=0.5) {
    C <- M/2
    if (opt.criterion=="OCI") {
        delta <- qnorm(1-alpha)+qnorm(beta)
        r1 <- RDEstimator(d, RDLFFunction(d, C, delta))
        r2 <- RDEstimator(d, RDLFFunction(d, C, 2*delta))
        return(r2$omega/(r1$delta*r1$sd+r1$omega))
    } else {
        ## From proof of Pratt result, it follows that the expected length is
        ## int pnorm(z_{1-alpha}-delta_t) dt, where delta_t is value of inverse
        ## two-class modulus omega{delta, F, 1}, or omega{delta, 1, F},
        ## depending on whether t > 0 or < 0.
        deltat <- function(t) sqrt(Q(d, RDgbC(d, t, C))) # delta_t
        integrand <- function(t)
            stats::pnorm(stats::qnorm(1-alpha)-sapply(t, deltat))

        ## By symmetry, half-length is given by value of integral over R_+. The
        ## integrand equals 1-alpha at zero, need upper cutoff
        upper <- 10
        while(integrand(upper)>1e-10) upper <- 2*upper
        return(stats::integrate(integrand, 1e-6, upper)$value /
                   RDTOpt.fit(d, 2*C, opt.criterion="FLCI", alpha=alpha)$hl)
    }
}
