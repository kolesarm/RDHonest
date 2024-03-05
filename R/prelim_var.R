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
        sigma2[k, ] <- Jk / (Jk+weights[k])*
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
PrelimVar <- function(d, se.initial="EHW") {
    ## Pilot bandwidth: either IK/ROT, or else Silverman (for actually computing
    ## IK) for uniform kernel making sure this results in enough distinct values
    ## on either side of threshold so we don't have perfect fit

    if (class(d)=="IP") {
        hmin <- max(sort(unique(abs(d$X)))[2], sort(abs(d$X))[4])
    } else {
        hmin <- max(sort(unique(d$X[d$p]))[3], sort(abs(unique(d$X[d$m])))[3],
                    sort(d$X[d$p])[4], sort(abs(d$X[d$m]))[4])
    }
    ## Use reduced form for FRD bandwidth selector
    drf <- d
    if (class(d)=="FRD") {
        drf$Y <- drf$Y[, 1]
        class(drf) <- "SRD"
    }

    if (se.initial == "EHW") {
        h1 <- if (class(d)=="IP") ROTBW(drf) else IKBW(drf)
        if (is.nan(h1)) {
            warning("Preliminary bandwidth is NaN, setting it to Inf")
            h1 <- Inf
        }
        r1 <- NPReg(d, max(h1, hmin), se.method="EHW")
    } else if (class(d) == "SRD" && se.initial == "Silverman") {
        ## Silverman only for RD/IK
        h1 <- max(1.84*stats::sd(d$X)/sum(length(d$X))^(1/5), hmin)
        r1 <- NPReg(d, h1, "uniform", order=0, se.method="EHW")
        ## Variance adjustment for backward compatibility
        lp <- length(r1$sigma2[d$p & r1$est_w != 0])
        lm <- length(r1$sigma2[d$m & r1$est_w != 0])
        r1$sigma2[d$p] <- r1$sigma2[d$p] * lp / (lp-1)
        r1$sigma2[d$m] <- r1$sigma2[d$m] * lm / (lm-1)
    } else {
        stop("This method for preliminary variance estimation not supported")
    }
    if (!is.null(d$clusterid))
        d$rho <- Moulton(r1$res, d)

    if (class(d) == "IP") {
        d$sigma2 <- rep(mean(r1$sigma2[r1$est_w != 0]), length(d$X))
    } else if (class(d)=="SRD") {
        d$sigma2 <- rep(NA, length(d$X))
        d$sigma2[d$p] <- mean(r1$sigma2[d$p & r1$est_w != 0])
        d$sigma2[d$m] <- mean(r1$sigma2[d$m & r1$est_w != 0])
    } else {
        d$sigma2 <- matrix(NA, nrow=length(d$X), ncol=4)
        d$sigma2[d$p, ] <- matrix(rep(colMeans(r1$sigma2[d$p & r1$est_w != 0, ]),
                                      each=sum(d$p)), nrow=sum(d$p))
        d$sigma2[d$m, ] <- matrix(rep(colMeans(r1$sigma2[d$m & r1$est_w != 0, ]),
                                      each=sum(d$m)), nrow=sum(d$m))
    }
    d
}

## Moulton estimate of rho, set rho=0 if no clustering
Moulton <- function(u, d) {
    m <- function(u, clusterid) {
        den <- sum(tapply(u[, 1], clusterid, length)^2)-NROW(u)
        if (den>0) {
            us <- apply(u, 2, function(x) tapply(x, clusterid, sum))
            as.vector(crossprod(us)-crossprod(u)) / den
        } else {
            rep(0, NCOL(u)^2)
        }
    }
    u <- as.matrix(u)
    ix <- !is.na(u[, 1])

    m(u[ix, , drop=FALSE], d$clusterid[ix])
}


## Rule of thumb bandwidth for inference at a point. Only used by PrelimVar
##
## Calculate bandwidth for inference at a point with local linear regression
## using method in Fan and Gijbels (1996, Chapter 4.2).
ROTBW <- function(d, kern="triangular") {
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
    s <- kernC[kernC$kernel==kern & kernC$order==order &
                   kernC$boundary==boundary, ]
    nu0 <- s$nu0
    mup <- s[[paste0("mu", order+1)]]

    ## STEP 3: Plug in
    B <- deriv * mup
    V <- sigma2 * nu0 /f0

    (V / (B^2 * 2 * (order+1) * N))^(1 / (2*order+3))
}

## Imbens and Kalyanaraman bandwidth. Only used by PrelimVar
##
##  Reproduce bandwidth from Section 6.2 in Imbens and Kalyanaraman (2012)
IKBW <- function(d, kern="triangular", verbose=FALSE) {
    Nm <- sum(d$m)
    Np <- sum(d$p)
    N <- Nm+Np

    ## STEP 0: Kernel constant
    s <- kernC[kernC$order==1 & kernC$boundary==TRUE & kernC$kernel==kern, ]
    const <- (s$nu0/s$mu2^2)^(1/5)

    ## STEP 1: Estimate f(0), sigma^2_(0) and sigma^2_+(0), using Silverman
    ## pilot bandwidth for uniform kernel
    d <- PrelimVar(d, se.initial="Silverman")
    h1 <- 1.84*stats::sd(d$X)/N^(1/5)
    f0 <- sum(abs(d$X) <= h1) / (2*N*h1)
    varm <- d$sigma2[d$m][1]
    varp <- d$sigma2[d$p][1]

    ## STEP 2: Estimate second derivatives m_{+}^(2) and m_{-}^(2)

    ## Estimate third derivative using 3rd order polynomial: Equation (14)
    m3 <- 6*stats::coef(stats::lm(d$Y ~ I(d$X>=0) + d$X + I(d$X^2) +
                                      I(d$X^3)))[5]

    ## Left and right bandwidths, Equation (15) and page 956.
    ## Optimal constant based on one-sided uniform Kernel, 7200^(1/7),
    h2m <- 7200^(1/7) * (varm / (f0*m3^2))^(1/7) * Nm^(-1/7)
    h2p <- 7200^(1/7) * (varp / (f0*m3^2))^(1/7) * Np^(-1/7)


    ## estimate second derivatives by local quadratic
    m2m <- 2*stats::coef(stats::lm(d$Y ~ d$X + I(d$X^2),
                                   subset = (d$X >= -h2m & d$X<0)))[3]
    m2p <- 2*stats::coef(stats::lm(d$Y ~ d$X + I(d$X^2),
                                   subset = (d$X <= h2p & d$X>=0)))[3]

    ## STEP 3: Calculate regularization terms, Equation (16)
    rm <- 2160*varm / (sum(d$X >= -h2m & d$X<0) * h2m^4)
    rp <- 2160*varp / (sum(d$X <= h2p & d$X>=0) * h2p^4)

    if (verbose)
        cat("\n h1: ", h1, "\n N_{-}, N_{+}: ", Nm, Np, "\n f(0): ", f0,
            "\n sigma^2_{+}(0): ", sqrt(varp), "^2\n sigma^2_{+}(0):",
            sqrt(varm), "^2", "\n m3: ", m3, "\n h_{2, +}:", h2p, "h_{2, -}:",
            h2m, "\n m^{(2)}_{+}: ", m2p, "m^{(2)}_{-}: ", m2m, "\n r_{+}:", rp,
            "\n r_{-}:", rm, "\n\n")

    ## Final bandwidth: Equation (17)
    unname(const * ((varp+varm) / (f0*N * ((m2p-m2m)^2+rm+rp)))^(1/5))
}
