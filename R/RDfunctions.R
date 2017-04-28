#' Local polynomial regression at a point, normalized to 0
#'
#' Calculate estimate of a function at a point and its variance given a
#' bandwidth using local polynomial regression of order \code{order}.
#'
#' @section Note: Nearest neighbor method assumes data are sorted so that
#'     \code{X[i] <= X[i+1]}
#' @param Y,X Outcome variable and regressor
#' @param h Bandwidth
#' @param K Kernel function
#' @param order Order of local regression 1 for linear, 2 for quadratic, etc.
#' @param sigma2 Optionally, supply estimates of \eqn{\sigma^{2}_{i}} (for
#'     \code{"supplied.var"} \code{se.method})
#' @template RDse
#' @keywords internal
LPReg <- function(X, Y, h, K, order=1, se.method=NULL, sigma2, J=3) {
    R <- outer(X, 0:order, "^")
    W <- K(X/h)
    Gamma <- crossprod(R, W*R)

    if (sum(W>0)<=order | h==0 | det(Gamma)==0)
        return(list(theta=0, sigma2=NA, var=NA, w=function(u) 0, eff.obs=0))

    beta <- drop(solve(Gamma, crossprod(R, W*Y)))
    hsigma2 <- drop(Y - R %*% beta)^2     # squared residuals

    EHW <- function(sigma2)
        (solve(Gamma, crossprod(sigma2*W*R, W*R)) %*% solve(Gamma))[1, 1]

    ## Robust variance-based formulae
    v <- c(NA, NA, NA, NA)
    names(v) <- c("supplied.var", "nn", "EHW", "demeaned")

    if("EHW" %in% se.method) v[3] <- EHW(hsigma2)
    if("supplied.var" %in% se.method) v[1] <- EHW(sigma2)
    if("demeaned" %in% se.method) {
        hsigma2 <- (Y-beta[1])^2
        v[4] <- EHW(hsigma2)
    }
    if("nn" %in% se.method) {
        hsigma2 <- sigmaNN(X, Y, J=J)
        v[2] <- EHW(hsigma2)
    }

    ## weight function if we think of the estimator as linear estimator
    w <- function(x) solve(Gamma, t(outer(x, 0:order, "^")*K(x/h)))[1, ]
    eff.obs <- 1/sum(w(X)^2)
    list(theta=beta[1], sigma2=hsigma2, var=v, w=w, eff.obs=eff.obs)
}


#' Local Polynomial Regression in Sharp RD
#'
#' Calculate sharp RD estimate and its variance given a bandwidth using local
#' polynomial regression of order \code{order}.
#'
#' @param d object of class \code{"RDData"}
#' @template RDBW
#' @template Kern
#' @template RDse
#' @param no.warning Don't warn about too few observations
#' @return list with elements:
#'
#' \describe{
#'     \item{estimate}{point estimate}
#'     \item{se}{Named vector of standard error estimates, as specified
#'               by \code{se.method}.}
#'     \item{wm,wp}{Implicit weight functions used}
#'
#'     \item{sigma2m, sigma2p}{Estimates of \eqn{sigma^2(X)} for values
#'             of running variable below (above) cutoff and receiving positive
#'             kernel weight. By default, estimates are based on squared
#'             regression residuals, as used in "EHW". If \code{"demeaned"}
#'             or \code{"nn"} is specifed,  estimates are based on that method,
#'             with \code{"nn"} method used if both are specified.}
#'     \item{eff.obs}{Number of effective observations}
#' }
#' @export
RDLPreg <- function(d, hp, kern="triangular", order=1, hm=hp, se.method="nn",
                    no.warning=FALSE, J=3) {

    K <- if (!is.function(kern)) {
             EqKern(kern, boundary=FALSE, order=0)
         } else {
             kern
         }
    Wm <- if (hm<=0) 0*d$Xm else K(d$Xm/hm) # kernel weights
    Wp <- if (hp<=0) 0*d$Xp else K(d$Xp/hp)

    ## variance calculations are faster if we only keep data with positive
    ## kernel weights
    d$Xp <- d$Xp[Wp>0]
    d$Xm <- d$Xm[Wm>0]

    if ((length(d$Xm) < 3*order | length(d$Xp) < 3*order)
        & no.warning==FALSE) {
        warning("Too few observations to compute RD estimates.\nOnly ",
                length(d$Xm), " control and ", length(d$Xp),
                " treated units with positive weights")
    }
    if ((length(unique(d$Xm)) <= order | length(unique(d$Xp)) <= order)
         & no.warning==FALSE) {
        warning("Too few distinct values to compute RD estimates.\nOnly ",
                length(unique(d$Xm)), " unique control and ",
                length(unique(d$Xp)), " unique treated values for ",
                "running variable with positive weights")
    }

    rm <- LPReg(d$Xm, d$Ym[Wm>0], hm, K, order, se.method, d$sigma2m[Wm>0], J)
    rp <- LPReg(d$Xp, d$Yp[Wp>0], hp, K, order, se.method, d$sigma2p[Wp>0], J)
    theta <- rp$theta-rm$theta

    ## Plug-in Asymptotic variance
    plugin <- NA
    if ("plugin" %in% se.method) {
        if (is.function(kern)) {
            ke <- EqKern(kern, boundary=TRUE, order=order)
            nu0 <- KernMoment(ke, moment=0, boundary=TRUE, "raw2")
        } else {
            nu0 <- RDHonest::kernC[RDHonest::kernC$kernel==kern &
                                   RDHonest::kernC$order==order &
                                   RDHonest::kernC$boundary==TRUE, "nu0"]
        }
        ## note we kept original outcomes, but only kept X's receiving positive
        ## weight
        N <- length(d$Yp) + length(d$Ym)
        f0 <- sum(length(d$Xp)+length(d$Xm)) / (N*(hp+hm))
        plugin <- nu0*(mean(rp$sigma2)/hp+mean(rm$sigma2)/hm) / (N*f0)
    }
    names(plugin) <- "plugin"
    se <- sqrt(c(rm$var+rp$var, plugin))
    ## eo <- 1/(1/rm$eff.obs+1/rp$eff.obs)
    eo <- rm$eff.obs+rp$eff.obs
    list(estimate=theta, se=se, wm=rm$w, wp=rp$w, sigma2m=rm$sigma2,
         sigma2p=rp$sigma2, eff.obs=eo)
}


#' method assumes X is sorted
#' @param J number of nearest neighbors
#' @keywords internal
sigmaNN <- function(X, Y, J=3) {
    n <- length(X)
    sigma2 <- vector(length=n)

    for (k in 1:n) {
        ## d is distance to Jth neighbor, exluding oneself
        s <- max(k-J, 1):k
        d <- sort(abs(c(X[s[-length(s)]], X[k:min(k+J, n)][-1])-X[k]))[J]
        ind <- (abs(X-X[k])<=d)
        ind[k] <- FALSE                 # exclude oneself
        Jk <- sum(ind)
        sigma2[k] <- Jk/(Jk+1)*(Y[k]-mean(Y[ind]))^2
    }
    sigma2
}

#' Compute preliminary estimate of variance
#'
#' @param d object of class \code{"RDData"}
#' @template RDseInitial
#' @return object of class \code{"RDData"} containing estimated variances.
#' @export
RDprelimVar <- function(d, se.initial="IKEHW") {
    ## Silverman pilot bandwidth for uniform kernel, making sure this
    ## results in enough distinct values on either side of threshold
    X <- c(d$Xm, d$Xp)
    h1 <- max(1.84*stats::sd(X)/sum(length(X))^(1/5),
              sort(unique(d$Xp))[2],
              abs(sort(unique(d$Xm), decreasing=TRUE)[2]))

    if (se.initial=="Silverman") {
        ## Equation (12) in IK
        d$sigma2m <- rep(stats::var(d$Ym[d$Xm >= -h1]), length(d$Xm))
        d$sigma2p <- rep(stats::var(d$Yp[d$Xp <= h1]), length(d$Xp))
    } else if (se.initial=="SilvermanNN") {
        ## This is how the Honest paper used to do it
        r1 <- RDLPreg(d=d, hp=h1, kern="uniform", order=1, se.method="nn")
        d$sigma2p <- rep(mean(r1$sigma2p), length(d$Xp))
        d$sigma2m <- rep(mean(r1$sigma2m), length(d$Xm))
    ## } else if (se.initial=="SilvermanEHW") {
    ##     r1 <- RDLPreg(d=d, hp=h1, kern="uniform", order=1, se.method="EHW")
    ##     d$sigma2p <- rep(mean(r1$sigma2p), length(d$Xp))
    ##     d$sigma2m <- rep(mean(r1$sigma2m), length(d$Xm))
    } else if (se.initial=="IKdemeaned") {
        r1 <- RDLPreg(d, IKBW.fit(d), se.method="demeaned")
        d$sigma2m <- rep(mean(r1$sigma2m), length(d$Xm))
        d$sigma2p <- rep(mean(r1$sigma2p), length(d$Xp))
    } else if (se.initial=="IKEHW"){
        r1 <- RDLPreg(d, IKBW.fit(d), se.method="EHW")
        d$sigma2m <- rep(mean(r1$sigma2m), length(d$Xm))
        d$sigma2p <- rep(mean(r1$sigma2p), length(d$Xp))
    } else {
        stop("Unknown method for estimating initial variance")
    }

    d
}
