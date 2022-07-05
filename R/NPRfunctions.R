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


## Local polynomial regression at a point, normalized to 0
## Calculate estimate of a function at a point and its variance given a
## bandwidth using local polynomial regression of order \code{order}.
## @section Note: Nearest neighbor method assumes data are sorted so that
##     \code{X[i] <= X[i+1]}
## @param Y,X Outcome variable and regressor
## @param h Bandwidth
## @param K Kernel function
## @param order Order of local regression 1 for linear, 2 for quadratic, etc.
## @param sigma2 Optionally, supply estimates of \eqn{\sigma^{2}_{i}} (for
##     \code{"supplied.var"} \code{se.method})
## @template RDse
LPReg <- function(X, Y, h, K, order=1, se.method=NULL, sigma2, J=3,
                  weights=rep(1L, length(X))) {
    R <- outer(X, 0:order, "^")
    W <- K(X/h)*weights
    Gamma <- crossprod(R, W*R)

    ## if (sum(W>0) <= order | h==0 |
    ##     class(try(solve(Gamma), silent=TRUE)) != "matrix")
    if (sum(W>0) <= order || h==0 ||
        inherits(try(solve(Gamma), silent=TRUE), "try-error"))
        return(list(theta=0, sigma2=NA, var=NA, w=0, eff.obs=0))

    ## weights if we think of the estimator as linear estimator
    ## w <- function(x) solve(Gamma, t(outer(x, 0:order, "^")*K(x/h)))[1, ]
    ## No longer works for weighted data
    w <- (W*R %*% solve(Gamma))[, 1]

    beta <- solve(Gamma, crossprod(R, W*Y))

    ## Squared residuals, allowing for multivariate Y
    sig <- function(r) {
        r[, rep(seq_len(ncol(r)), each=ncol(r))] *
        r[, rep(seq_len(ncol(r)), ncol(r))]
    }
    hsigma2 <- sig(Y - R %*% beta)

    ## Robust variance-based formulae
    v <- matrix(nrow=NCOL(hsigma2), ncol=3)
    colnames(v) <- c("supplied.var", "nn", "EHW")
    EHW <- function(sigma2) {
        if (NCOL(sigma2)==1)
            sum(w^2 * sigma2)
        else
            colSums(w^2 * sigma2)
    }
    if("EHW" %in% se.method) v[, 3] <- EHW(hsigma2)
    if("supplied.var" %in% se.method) v[, 1] <- EHW(sigma2)
    if("nn" %in% se.method) {
        hsigma2 <- sigmaNN(X, Y, J=J, weights)
        v[, 2] <- EHW(hsigma2)
    }

    list(theta=beta[1, ], sigma2=hsigma2, var=drop(v),
         w=w, eff.obs=1/sum(w^2))
}


## Nonparametric Regression
##
## Calculate fuzzy or sharp RD estimate, or estimate of a conditional mean at a
## point (depending on the class of \code{d}), and its variance using local
## polynomial regression of order \code{order}.
##
## @param d object of class \code{"LPPData"}, \code{"RDData"}, or
##     \code{"FRDData"}
## @template RDBW
## @template Kern
## @template RDse
## @param no.warning Don't warn about too few observations
## @return list with elements:
##
## \describe{
##     \item{estimate}{point estimate}
##     \item{se}{Named vector of standard error estimates, as specified
##               by \code{se.method}.}
##     \item{w}{Implicit weight function used}
##
##     \item{sigma2}{Estimate of \eqn{\sigma^2(X)}{sigma^2(X)} for values of
##             \eqn{X} receiving positive kernel weight. By default, estimates
##             are based on squared regression residuals, as used in
##             \code{"EHW"}. If \code{se.method="nn"} is specified, estimates
##             are based on that method, with \code{"nn"} method used if both
##             are specified.}
##
##      \item{eff.obs}{Number of effective observations}
##
## }
NPRreg.fit <- function(d, h, kern="triangular", order=1, se.method="nn",
                       no.warning=FALSE, J=3) {

    K <- kern
    if (!is.function(kern))
        K <- EqKern(kern, boundary=FALSE, order=0)

    if (inherits(d, "LPPData")) {
        ## Keep only positive kernel weights
        W <- if (h<=0) 0*d$X else K(d$X/h) # kernel weights
        d$X <- d$X[W>0]
    } else {
        Wm <- if (h<=0) 0*d$Xm else K(d$Xm/h)
        Wp <- if (h<=0) 0*d$Xp else K(d$Xp/h)
        d$Xp <- d$Xp[Wp>0]
        d$Xm <- d$Xm[Wm>0]
        if(!is.null(d$sigma2p)) d$sigma2p <- as.matrix(d$sigma2p)
        if(!is.null(d$sigma2m)) d$sigma2m <- as.matrix(d$sigma2m)
    }
    nob <- max(length(d$X), min(length(d$Xp), length(d$Xm)))
    uob <- max(length(unique(d$X)),
               min(length(unique(d$Xp)), length(unique(d$Xm))))

    if (no.warning==FALSE && (nob <= 3*order || uob <= order))
        warning("Too few observations to compute estimates.\nOnly ",
                nob, " units with positive weights and ",
                uob, " unique values for ",
                "independent variable with positive weights")

    if (inherits(d, "LPPData")) {
        r <- LPReg(d$X, d$Y[W>0], h, K, order, se.method, d$sigma2[W>0],
                   J, weights=d$w[W>0])
        ## Estimation weights
        W[W>0] <- r$w
        return(list(estimate=r$theta, se=sqrt(r$var), w=W,
                    sigma2=r$sigma2, eff.obs=r$eff.obs))
    }
    rm <- LPReg(d$Xm, as.matrix(d$Ym)[Wm>0, ], h, K, order, se.method,
                d$sigma2m[Wm>0, ], J, weights=d$wm[Wm>0])
    rp <- LPReg(d$Xp, as.matrix(d$Yp)[Wp>0, ], h, K, order, se.method,
                d$sigma2p[Wp>0, ], J, weights=d$wp[Wp>0])
    Wm[Wm>0] <- rm$w
    Wp[Wp>0] <- rp$w
    ret <- list(estimate=NULL, se=NULL, wm=Wm, wp=Wp, sigma2m=rm$sigma2,
                sigma2p=rp$sigma2, eff.obs=rm$eff.obs+rp$eff.obs)
    if (inherits(d, "RDData")) {
        ret$estimate <- rp$theta-rm$theta
        ret$se <- sqrt(rm$var+rp$var)
    } else if (inherits(d, "FRDData")) {
        ret$fs <- rp$theta[2]-rm$theta[2]
        ret$estimate <- (rp$theta[1]-rm$theta[1]) / ret$fs
        ret$se <- sqrt(drop(t(rp$var+rm$var) %*%
                            c(1, -ret$estimate, -ret$estimate, ret$estimate^2)
                            ) / ret$fs^2)

    }
    ret
}

## Compute preliminary estimate of variance: nn, EHW, Silverman
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
    } else if (se.initial == "Silverman") {
        X <- if (inherits(d, "LPPData")) d$X else c(d$Xm, d$Xp)
        Xmin <- if (inherits(d, "LPPData")) {
                    sort(abs(unique(d$X)))[2]
                } else  {
                    max(sort(unique(d$Xp))[2], sort(abs(unique(d$Xm)))[2])
                }
        h1 <- max(1.84*stats::sd(X)/sum(length(X))^(1/5), Xmin)
        r1 <- NPRreg.fit(d=d, h=h1, kern="uniform", order=0, se.method="EHW")
        ## Variance adjustment for backward compatibility
        if (inherits(d, "LPPData")) {
            r1$sigma2 <- r1$sigma2*length(r1$sigma2) / (length(r1$sigma2)-1)
        } else {
            r1$sigma2p <- r1$sigma2p*length(r1$sigma2p) /
                (length(r1$sigma2p)-1)
            r1$sigma2m <- r1$sigma2m*length(r1$sigma2m) /
                (length(r1$sigma2m)-1)
        }
    } else {
        stop("Unknown method for estimating initial variance")
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

## Rule of thumb for choosing M
##
## Use global quartic regression to estimate a bound on the second derivative
## for inference under under second order Hölder class. For RD, use a separate
## regression on either side of the cutoff
NPR_MROT.fit <- function(d) {
    if (inherits(d, "RDData")) {
        max(NPR_MROT.fit(LPPData(data.frame(Y=d$Yp, X=d$Xp), 0)),
            NPR_MROT.fit(LPPData(data.frame(Y=d$Ym, X=d$Xm), 0)))
    } else if (inherits(d, "FRDData")) {
        c(M1=max(NPR_MROT.fit(LPPData(data.frame(Y=d$Yp[, 1], X=d$Xp), 0)),
                 NPR_MROT.fit(LPPData(data.frame(Y=d$Ym[, 1], X=d$Xm), 0))),
          M2=max(NPR_MROT.fit(LPPData(data.frame(Y=d$Yp[, 2], X=d$Xp), 0)),
                 NPR_MROT.fit(LPPData(data.frame(Y=d$Ym[, 2], X=d$Xm), 0))))
    } else if (inherits(d, "LPPData")) {
        ## STEP 1: Estimate global polynomial regression
        r1 <- unname(stats::lm(d$Y ~ 0 + outer(d$X, 0:4, "^"))$coefficients)
        f2 <- function(x) abs(2*r1[3]+6*x*r1[4]+12*x^2*r1[5])
        ## maximum occurs either at endpoints, or else at the extremum,
        ## -r1[4]/(4*r1[5]), if the extremum is in the support
        f2e <- if(abs(r1[5])<=1e-10) Inf else -r1[4]/(4*r1[5])
        M <- max(f2(min(d$X)), f2(max(d$X)))
        if (min(d$X) < f2e && max(d$X) > f2e) M <- max(f2(f2e), M)

        M
    }
}
