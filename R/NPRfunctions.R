## Nonparametric Regression
##
## Calculate fuzzy or sharp RD estimate, or estimate of a conditional mean at a
## point (depending on the class of \code{d}), and its variance using local
## polynomial regression of order \code{order}.
NPReg <- function(d, h, kern="triangular", order=1, se.method="nn", J=3) {
    if (!is.function(kern))
        kern <- EqKern(kern, boundary=FALSE, order=0)
    ## Keep only positive kernel weights
    X <- drop(d$X)
    W <- if (h<=0) 0*X else kern(X/h)*d$w # kernel weights
    ## if (!is.null(d$sigma2)) d$sigma2 <- as.matrix(d$sigma2)
    Z <- outer(X, 0:order, "^")
    nZ <- c("(Intercept)", colnames(d$X), paste0("I(", colnames(d$X), "^",
                                                 seq_len(order), ")")[-1])
    colnames(Z) <- nZ[seq_len(order+1)]
    if (!inherits(d, "IP")) {
        ZZ <- (X>=0)*Z
        colnames(ZZ) <- c(paste0("I(", colnames(d$X), ">0)"),
                          paste0(paste0("I(", colnames(d$X), ">0):"),
                                 colnames(Z))[-1])
        Z <- cbind(ZZ, Z, d$covs)
    }
    r0 <- stats::lm.wfit(x=Z, y=d$Y, w=W)
    class(r0) <- c(if (ncol(d$Y)>1) "mlm", "lm")
    if (any(is.na(r0$coefficients))) {
        return(list(estimate=0, se=NA, est_w=W*0, sigma2=NA*d$Y, eff.obs=0,
                    fs=NA, lm=r0))
    }
    wgt <- W*0
    ok <- W!=0
    wgt[ok] <- solve(qr.R(r0$qr), t(sqrt(W[W>0])*qr.Q(r0$qr)))[1, ]
    ## To compute effective observations, rescale against uniform kernel
    Wu <- d$w * (abs(X)<=h)
    q_u <- qr(sqrt(Wu) * Z)
    wgt_u <- solve(qr.R(q_u), t(sqrt(Wu) * qr.Q(q_u)))[1, ]
    eff.obs <- sum(Wu)*sum(wgt_u^2/d$w)/sum(wgt^2/d$w)

    ny <- NCOL(r0$residuals)
    ## Squared residuals, allowing for multivariate Y
    HC <- function(r) r[, rep(seq_len(ny), each=ny)] * r[, rep(seq_len(ny), ny)]

    ## For RD, compute variance separately on either side of cutoff
    ## TODO: better implement
    NN <- function(X) {
        res <- matrix(0, nrow=length(X), ncol=ny^2)
        res[ok] <-
            if (!inherits(d, "IP"))
                rbind(as.matrix(sigmaNN(X[d$m & ok], d$Y[d$m & ok, ], J,
                                        d$w[d$m & ok])),
                      as.matrix(sigmaNN(X[d$p & ok], d$Y[d$p & ok, ], J,
                                        d$w[d$p & ok])))
            else
                sigmaNN(X[ok], d$Y[ok, ], J, d$w[ok])
        res
    }

    hsigma2 <- switch(se.method, nn=NN(X), EHW=HC(as.matrix(r0$residuals)),
                      supplied.var=as.matrix(d$sigma2))

    ## Variance
    if (is.null(d$clusterid)) {
        V <- colSums(as.matrix(wgt^2 * hsigma2))
    } else if (se.method == "supplied.var") {
        V <- colSums(as.matrix(wgt^2 * hsigma2))+
            d$rho * (sum(tapply(wgt, d$clusterid, sum)^2)-sum(wgt^2))
    } else {
        us <- apply(as.matrix(wgt*r0$residuals)[ok, , drop=FALSE], 2,
                    function(x) tapply(x, d$clusterid[ok], sum))
        V <- as.vector(crossprod(us))
    }
    ret <- list(estimate=r0$coefficients[1], se=sqrt(V[1]), est_w=wgt,
                sigma2=hsigma2, eff.obs=eff.obs, fs=NA, lm=r0)

    if (inherits(d, "FRD")) {
        ret$fs <- r0$coefficients[1, 2]
        ret$estimate <- r0$coefficients[1, 1]/r0$coefficients[1, 2]
        names(ret$estimate) <- colnames(r0$coefficients)[2]
        ret$se <- sqrt(sum(c(1, -ret$estimate, -ret$estimate,
                             ret$estimate^2) * V) / ret$fs^2)
    }
    ret
}


## Rule of thumb for choosing M
##
## Use global quartic regression to estimate a bound on the second derivative
## for inference under under second order Hölder class. For RD, use a separate
## regression on either side of the cutoff
MROT <- function(d) {
    if (inherits(d, "SRD")) {
        max(MROT(data.frame(Y=d$Y[d$p], X=d$X[d$p], w=d$w[d$p])),
            MROT(data.frame(Y=d$Y[d$m], X=d$X[d$m], w=d$w[d$m])))
    } else if (inherits(d, "FRD")) {
        c(MY=max(MROT(data.frame(Y=d$Y[d$p, 1], X=d$X[d$p], w=d$w[d$p])),
                 MROT(data.frame(Y=d$Y[d$m, 1], X=d$X[d$m], w=d$w[d$m]))),
          MD=max(MROT(data.frame(Y=d$Y[d$p, 2], X=d$X[d$p], w=d$w[d$p])),
                 MROT(data.frame(Y=d$Y[d$m, 2], X=d$X[d$m], w=d$w[d$m]))))
    } else {
        ## STEP 1: Estimate global polynomial regression
        r1 <- unname(stats::lm.wfit(y=d$Y, x=outer(drop(d$X), 0:4, "^"),
                                    w=d$w)$coefficients)
        if (length(unique(d$X))<4 || any(is.na(r1)))
            stop(paste0("Insufficient unique values of the running",
                        " variable to compute rule of thumb for M."))
        f2 <- function(x) abs(2*r1[3]+6*x*r1[4]+12*x^2*r1[5])
        ## maximum occurs either at endpoints, or else at the extremum,
        ## -r1[4]/(4*r1[5]), if the extremum is in the support
        f2e <- if (abs(r1[5])<=1e-10) Inf else -r1[4] / (4*r1[5])
        M <- max(f2(min(d$X)), f2(max(d$X)))
        if (min(d$X) < f2e && max(d$X) > f2e) M <- max(f2(f2e), M)

        M
    }
}
