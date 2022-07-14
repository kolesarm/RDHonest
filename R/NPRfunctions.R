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
LPReg <- function(X, Y, h, K, order=1, se.method=NULL, sigma2, J=3,
                  weights=rep(1L, length(X))) {
    R <- outer(X, 0:order, "^")
    W <- K(X/h)*weights
    Gamma <- crossprod(R, W*R)

    if (sum(W>0) <= order || h==0 ||
        inherits(try(solve(Gamma), silent=TRUE), "try-error"))
        return(list(theta=0, sigma2=NA, var=NA, w=0, eff.obs=0))

    ## weights if we think of the estimator as linear estimator
    wgt <- (W*R %*% solve(Gamma))[, 1]
    beta <- solve(Gamma, crossprod(R, W*Y))

    ## Squared residuals, allowing for multivariate Y
    sig <- function(r) {
        r[, rep(seq_len(ncol(r)), each=ncol(r))] *
        r[, rep(seq_len(ncol(r)), ncol(r))]
    }
    hsigma2 <- switch(se.method,
                      nn=sigmaNN(X, Y, J, weights),
                      EHW=sig(Y - R %*% beta),
                      supplied.var=sigma2)
    ## Robust variance-based formulae
    EHW <- function(sigma2) {
        if (NCOL(sigma2)==1)
            sum(wgt^2 * sigma2)
        else
            colSums(wgt^2 * sigma2)
    }

    ## eff.obs=1/sum(w^2) TODO
    list(theta=beta[1, ], sigma2=hsigma2, var=EHW(hsigma2), w=wgt, eff.obs=NA)
}


## Nonparametric Regression
##
## Calculate fuzzy or sharp RD estimate, or estimate of a conditional mean at a
## point (depending on the class of \code{d}), and its variance using local
## polynomial regression of order \code{order}.
NPRreg.fit <- function(d, h, kern="triangular", order=1, se.method="nn", J=3) {
    if (!is.function(kern))
        kern <- EqKern(kern, boundary=FALSE, order=0)
    ## Keep only positive kernel weights
    W <- if (h<=0) 0*d$X else kern(d$X/h) # kernel weights
    sigma2 <- matrix(NA, nrow=length(W), ncol=NCOL(d$Y)^2)
    if (d$class=="IP") {
        r <- LPReg(d$X[W>0], d$Y[W>0], h, kern, order, se.method, d$sigma2[W>0],
                   J, weights=d$w[W>0])
        ## Estimation weights
        sigma2[W>0] <- r$sigma2
        W[W>0] <- r$w
        return(list(estimate=r$theta, se=sqrt(r$var), w=W,
                    sigma2=sigma2,
                    eff.obs=r$eff.obs, fs=NA))
    }
    if(!is.null(d$sigma2)) d$sigma2 <- as.matrix(d$sigma2)
    rm <- LPReg(d$X[W>0 & d$m], as.matrix(d$Y)[W>0 & d$m, ], h, kern, order,
                se.method, d$sigma2[W>0 & d$m, ], J, weights=d$w[W>0 & d$m])
    rp <- LPReg(d$X[W> 0 & d$p], as.matrix(d$Y)[W>0 & d$p, ], h, kern, order,
                se.method, d$sigma2[W>0 & d$p, ], J, weights=d$w[W>0 & d$p])

    sigma2[W>0 & d$m, ] <- rm$sigma2
    sigma2[W>0 & d$p, ] <- rp$sigma2
    W[W>0 & d$m] <- rm$w
    W[W>0 & d$p] <- rp$w

    ret <- list(estimate=rp$theta[1]-rm$theta[1], se=sqrt(rm$var[1]+rp$var[1]),
                w=W, sigma2=drop(sigma2), eff.obs=rm$eff.obs+rp$eff.obs, fs=NA)

    if (d$class=="FRD") {
        ret$fs <- rp$theta[2]-rm$theta[2]
        ret$estimate <- (rp$theta[1]-rm$theta[1]) / ret$fs
        ret$se <- sqrt(sum(c(1, -ret$estimate, -ret$estimate, ret$estimate^2) *
                           (rp$var+rm$var)) / ret$fs^2)
    }
    ret
}


## Rule of thumb for choosing M
##
## Use global quartic regression to estimate a bound on the second derivative
## for inference under under second order Hölder class. For RD, use a separate
## regression on either side of the cutoff
NPR_MROT.fit <- function(d) {
    if (d$class=="SRD") {
        max(NPR_MROT.fit(NPRData(data.frame(Y=d$Y[d$p], X=d$X[d$p]), 0, "IP")),
            NPR_MROT.fit(NPRData(data.frame(Y=d$Y[d$m], X=d$X[d$m]), 0, "IP")))
    } else if (d$class=="FRD") {
        c(M1=max(NPR_MROT.fit(NPRData(data.frame(Y=d$Y[d$p, 1], X=d$X[d$p]), 0,
                                      "IP")),
                 NPR_MROT.fit(NPRData(data.frame(Y=d$Y[d$m, 1], X=d$X[d$m]), 0,
                                      "IP"))),
          M2=max(NPR_MROT.fit(NPRData(data.frame(Y=d$Y[d$p, 2], X=d$X[d$p]), 0,
                                      "IP")),
                 NPR_MROT.fit(NPRData(data.frame(Y=d$Y[d$m, 2], X=d$X[d$m]), 0,
                                      "IP"))))
    } else if (d$class=="IP") {
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
