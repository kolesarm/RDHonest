## Functions for variance estimation, including IK and ROT bandwidth selectors

## method assumes X is sorted
## @param J number of nearest neighbors
sigmaNN <- function(X, Y, J=3, weights=rep(1L, length(X))) {
    n <- length(X)
    sigma2 <- matrix(nrow=n, ncol=NCOL(Y)^2)

    for (k in seq_along(X)) {
        ## d is distance to Jth neighbor, exluding oneself
        s <- max(k-J, 1):k
        d <- sort(abs(c(X[s[-length(s)]], X[k:min(k+J, n)][-1])-X[k]))[J]
        ind <- (abs(X-X[k])<=d)
        ind[k] <- FALSE                 # exclude oneself
        Jk <- sum(weights[ind])
        sigma2[k, ] <- Jk/(Jk+weights[k])*
            if (NCOL(Y)>1)
                as.vector(outer(Y[k, ]-colSums(weights[ind]*Y[ind, ])/Jk,
                                Y[k, ]-colSums(weights[ind]*Y[ind, ])/Jk))
            else
                (Y[k]-sum(weights[ind]*Y[ind])/Jk)^2
    }
    drop(sigma2)
}


## Compute preliminary variance estimate: EHW or Silverman. Silverman used by IK
##
## Compute estimate of variance, which can then be used in optimal bandwidth
## calculations. These estimates are unweighted.
NPRPrelimVar.fit <- function(d, se.initial="EHW") {
    ## Pilot bandwidth: either IK/ROT, or else Silverman for uniform kernel,
    ## making sure this results in enough distinct values on either side of
    ## threshold
    if (se.initial == "EHW") {
        drf <- d
        ## Use reduced form for FRD
        if (inherits(d, "FRDData")) {
            drf$Yp <- drf$Yp[, 1]
            drf$Ym <- drf$Ym[, 1]
            class(drf) <- "RDData"
        }
        h1 <- if (inherits(d, "LPPData")) ROTBW.fit(drf) else IKBW.fit(drf)
        if(is.nan(h1)) {
            warning("Preliminary bandwidth is NaN, setting it to Inf")
            h1 <- Inf
        }
        r1 <- NPRreg.fit(d, h1, se.method=se.initial)
    } else if (!inherits(d, "LPPData")) { # Silverman only for RD/IK
        X <- c(d$Xm, d$Xp)
        Xmin <- max(sort(unique(d$Xp))[2], sort(abs(unique(d$Xm)))[2])
        h1 <- max(1.84*stats::sd(X)/sum(length(X))^(1/5), Xmin)
        r1 <- NPRreg.fit(d=d, h=h1, kern="uniform", order=0, se.method="EHW")
        ## Variance adjustment for backward compatibility
        r1$sigma2p <- r1$sigma2p*length(r1$sigma2p) / (length(r1$sigma2p)-1)
        r1$sigma2m <- r1$sigma2m*length(r1$sigma2m) / (length(r1$sigma2m)-1)
    }

    if (inherits(d, "LPPData")) {
        d$sigma2 <- rep(mean(r1$sigma2), length(d$X))
    } else if (inherits(d, "RDData")) {
        d$sigma2p <- rep(mean(r1$sigma2p), length(d$Xp))
        d$sigma2m <- rep(mean(r1$sigma2m), length(d$Xm))
    } else {
        d$sigma2m <- matrix(rep(colMeans(r1$sigma2m), each=length(d$Xm)),
                                nrow=length(d$Xm))
        d$sigma2p <- matrix(rep(colMeans(r1$sigma2p), each=length(d$Xp)),
                            nrow=length(d$Xp))
    }
    d
}


## Rule of thumb bandwidth for inference at a point. Only used by
## NPRPrelimVar.fit
##
## Calculate bandwidth for inference at a point with local linear regression
## using method in Fan and Gijbels (1996, Chapter 4.2).
ROTBW.fit <- function(d, kern="triangular") {
    X <- d$X
    boundary <- if ((min(X)>=0) || (max(X)<=0)) TRUE else FALSE
    N <- length(d$X)

    ## STEP 0: Estimate f_X(0) using Silverman
    h1 <- 1.843 *
        min(stats::sd(X), (stats::quantile(X, 0.75) -
                           stats::quantile(X, 0.25)) / 1.349) / N^(1/5)
    f0 <- sum(abs(X) <= h1) / (2*N*h1)

    ## STEP 1: Estimate (p+1)th derivative and sigma^2 using global polynomial
    ## regression
    order <- 1
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


## Imbens and Kalyanaraman bandwidth. Only used by NPRPrelimVar.fit
##
##  Reproduce bandwidth from Section 6.2 in Imbens and Kalyanaraman (2012)
IKBW.fit <- function(d, kern="triangular", verbose=FALSE) {
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
