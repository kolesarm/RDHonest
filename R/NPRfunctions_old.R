## Local polynomial regression/RD at a point, normalized to 0
## Only used by NPReg
## @param sigma2 Optionally, supply estimates of \eqn{\sigma^{2}_{i}} (for
##     \code{"supplied.var"} \code{se.method})
LPReg <- function(X, Y, h, K, order=1, se.method=NULL, sigma2, J=3,
                  weights=rep(1L, length(X)), RD=FALSE, rho=NULL,
                  clusterid=NULL, covs=matrix(nrow=NROW(X), ncol=0)) {
    X <- drop(X)
    R <- outer(X, 0:order, "^")
    if (RD)
        R <- cbind((X>=0)*R, R, covs)
    W <- K(X/h)*weights
    Gamma <- crossprod(R, W*R)

    if (h==0 || inherits(try(solve(Gamma), silent=TRUE), "try-error"))
        return(list(theta=0, sigma2=NA*Y, res=NA*Y,
                    var=NA*stats::var(Y), est_w=W*0, eff.obs=0))
    ## weights if we think of the estimator as linear estimator
    wgt <- (W*R %*% solve(Gamma))[, 1]
    beta <- solve(Gamma, crossprod(R, W*Y))

    ## Squared residuals, allowing for multivariate Y
    sig <- function(r) {
        r[, rep(seq_len(ncol(r)), each=ncol(r))] * r[, rep(seq_len(ncol(r)),
                                                           ncol(r))]
    }
    ## For RD, compute variance separately on either side of cutoff
    signn <- function(X) {
        if (RD && NCOL(Y)==1) # Sharp RD
            c(sigmaNN(X[X<0], Y[X<0], J, weights[X<0]),
              sigmaNN(X[X>=0], Y[X>=0], J, weights[X>=0]))
        else if (RD && NCOL(Y)>1) # Fuzzy RD
            rbind(sigmaNN(X[X<0], Y[X<0, ], J, weights[X<0]),
                  sigmaNN(X[X>=0], Y[X>=0, ], J, weights[X>=0]))
        else
            sigmaNN(X, Y, J, weights)
    }
    res <- Y - R %*% beta
    hsigma2 <- switch(se.method, nn=signn(X), EHW=sig(res), supplied.var=sigma2)

    ## To compute effective observations, rescale against uniform kernel
    wgt_unif <- (weights*R %*% solve(crossprod(R, weights*R)))[, 1]
    eff.obs <- sum(weights)*sum(wgt_unif^2/weights)/sum(wgt^2/weights)
    ## Variance
    if (is.null(clusterid)) {
        V <- colSums(as.matrix(wgt^2 * hsigma2))
    } else if (se.method == "supplied.var") {
        V <- colSums(as.matrix(wgt^2 * hsigma2))+
            rho * (sum(tapply(wgt, clusterid, sum)^2)-sum(wgt^2))
    } else {
        us <- apply(wgt*res, 2, function(x) tapply(x, clusterid, sum))
        V <- as.vector(crossprod(us))
    }
    list(theta=beta[1, ], sigma2=hsigma2, res=res, var=V, est_w=wgt,
         eff.obs=eff.obs,
         gamma=beta[seq_len(NCOL(covs))+NROW(beta)-NCOL(covs), ])
}




## Nonparametric Regression
##
## Calculate fuzzy or sharp RD estimate, or estimate of a conditional mean at a
## point (depending on the class of \code{d}), and its variance using local
## polynomial regression of order \code{order}.
NPReg2 <- function(d, h, kern="triangular", order=1, se.method="nn", J=3) {
    if (!is.function(kern))
        kern <- EqKern(kern, boundary=FALSE, order=0)
    ## Keep only positive kernel weights
    W <- if (h<=0) 0*d$X else kern(d$X/h)*d$w # kernel weights
    if (!is.null(d$sigma2)) d$sigma2 <- as.matrix(d$sigma2)
    r <- LPReg(d$X[W>0], as.matrix(d$Y)[W>0, ], h, kern, order, se.method,
               d$sigma2[W>0, ], J, weights=d$w[W>0], RD = !inherits(d, "IP"),
               d$rho, d$clusterid[W>0])
    sigma2 <- matrix(NA, nrow=length(W), ncol=NCOL(d$Y)^2)
    sigma2[W>0, ] <- r$sigma2
    res <- matrix(NA, nrow=length(W), ncol=NCOL(d$Y))
    res[W>0, ] <- r$res # residuals
    W[W>0] <- r$est_w

    ret <- list(estimate=r$theta[1], se=sqrt(r$var[1]),
                est_w=W, res=drop(res), w=d$w,
                sigma2=drop(sigma2), eff.obs=r$eff.obs, fs=NA)

    if (inherits(d, "FRD")) {
        ret$fs <- r$theta[2]
        ret$estimate <- r$theta[1]/r$theta[2]
        ret$se <- sqrt(sum(c(1, -ret$estimate, -ret$estimate,
                             ret$estimate^2) * r$var) / ret$fs^2)
    }
    ret
}
