#' Lower bound on smoothness constant M in sharp RD designs
#'
#' Estimate a lower bound on the smoothness constant M and provide a lower
#' confidence interval for it, using method described in supplement to
#' Kolesár and Rothe (2018).
#'
#' @param object An object of class \code{"RDResults"}, typically a result of a
#'     call to \code{\link{RDHonest}}.
#' @param s Number of support points that curvature estimates should average
#'     over.
#' @param sclass Smoothness class, either \code{"T"} for Taylor or \code{"H"}
#'     for Hölder class.
#' @param alpha determines confidence level \code{1-alpha}.
#' @param separate If \code{TRUE}, report estimates separately for data above
#'     and below cutoff. If \code{FALSE}, report pooled estimates.
#' @param multiple If \code{TRUE}, use multiple curvature estimates. If
#'     \code{FALSE}, only use a single curvature estimate using observations
#'     closest to the cutoff.
#' @return Returns a data frame wit the following columns:
#'
#' \describe{

#' \item{\code{estimate}}{Point estimate for lower bounds for M. }
#'
#' \item{\code{conf.low}}{Lower endpoint for a one-sided confidence interval
#'                        for M}
#'
#' }
#'
#' The data frame has a single row if \code{separate==FALSE}; otherwise it has
#' two rows, corresponding to smoothness bound estimates and confidence
#' intervals below and above the cutoff, respectively.
#' @references{
#'
#' \cite{Michal Kolesár and Christoph Rothe. Inference in regression
#'       discontinuity designs with a discrete running variable.
#'       American Economic Review, 108(8):2277—-2304,
#'       August 2018. \doi{10.1257/aer.20160945}}
#'
#' }
#' @examples
#' ## Subset data to increase speed
#' r <- RDHonest(log(earnings)~yearat14, data=cghs,
#'               subset=abs(yearat14-1947)<10,
#'               cutoff=1947, M=0.04, h=3)
#' RDSmoothnessBound(r, s=2)
#' @export
RDSmoothnessBound <- function(object, s, separate=FALSE, multiple=TRUE,
                              alpha=0.05, sclass="H") {
    d <- PrelimVar(object$data, se.initial="EHW")

    ## Curvature estimate based on jth set of three points closest to zero
    Dk <- function(Y, X, xu, s2, j) {
        I1 <- X >= xu[3*j*s-3*s+1]  & X <= xu[3*j*s-2*s]
        I2 <- X >= xu[3*j*s-2*s+1]  & X <= xu[3*j*s-s]
        I3 <- X >= xu[3*j*s-  s+1]  & X <= xu[3*j*s]

        lam <- (mean(X[I3])-mean(X[I2])) / (mean(X[I3])-mean(X[I1]))
        if  (sclass=="T") {
            den <- (1-lam)*mean(X[I3]^2) + lam*mean(X[I1]^2) + mean(X[I2]^2)
        } else {
            den <- (1-lam)*mean(X[I3]^2) + lam*mean(X[I1]^2) - mean(X[I2]^2)
        }
        ## Delta is lower bound on M by lemma S2 in Kolesar and Rothe
        Del <- 2 * (lam*mean(Y[I1]) + (1-lam) * mean(Y[I3])-mean(Y[I2])) / den
        ## Variance of Delta
        VD <- 4 * (lam^2*mean(s2[I1])/sum(I1) + (1-lam)^2*mean(s2[I3]) /
                       sum(I3) + mean(s2[I2])/sum(I2)) / den^2
        c(Del, sqrt(VD), mean(Y[I1]), mean(Y[I2]), mean(Y[I3]), range(X[I1]),
          range(X[I2]), range(X[I3]))
    }

    xp <- unique(d$X[d$p])
    xm <- sort(unique(abs(d$X[d$m])))

    Dpj <- function(j) Dk(d$Y[d$p], d$X[d$p], xp, d$sigma2[d$p], j)
    Dmj <- function(j) Dk(d$Y[d$m], abs(d$X[d$m]), xm, d$sigma2[d$m], j)

    Sp <- floor(length(xp) / (3*s))
    Sm <- floor(length(xm) / (3*s))
    if (min(Sp, Sm) == 0)
        stop("Value of s is too big")

    if (!multiple) {
        Sp <- Sm <- 1
    }
    Dp <- vapply(seq_len(Sp), Dpj, numeric(11))
    Dm <- vapply(seq_len(Sm), Dmj, numeric(11))

    ## Critical value
    cv <- function(M, Z, sd, alpha) {
        if (ncol(Z)==1) {
            CVb(M/sd, alpha=alpha)
        } else {
            S <- Z+M*outer(rep(1, nrow(Z)), 1/sd)
            maxS <- abs(S[cbind(seq_len(nrow(Z)), max.col(S))])
            unname(stats::quantile(maxS, 1-alpha))
        }
    }

    hatM <- function(D) {
        ts <- abs(D[1, ]/D[2, ]) # sup_t statistic
        maxt <- D[, which.max(ts)]

        Z <- matrix(stats::rnorm(10000*ncol(D)), ncol=ncol(D))
        ## median unbiased point estimate and lower CI
        hatm <- lower <- 0
        if (max(ts) > cv(0, Z, D[2, ], 1/2)) {
            hatm <- FindZero(function(m) max(ts)-cv(m, Z, D[2, ], 1/2),
                             negative=FALSE)
        }
        if (max(ts) >= cv(0, Z, D[2, ], alpha)) {
            lower <- FindZero(function(m) max(ts)-cv(m, Z, D[2, ], alpha),
                              negative=FALSE)
        }
        list(estimate=hatm, conf.low=lower,
             diagnostics=c(Delta=maxt[1], sdDelta=maxt[2],
                           y1=maxt[3], y2=maxt[4], y3=maxt[5],
                           I1=maxt[6:7], I2=maxt[8:9], I3=maxt[10:11]))
    }

    if (separate) {
        withr::with_seed(42, po <- hatM(Dp))
        withr::with_seed(42, ne <- hatM(Dm))
        ret <- data.frame(rbind("Below cutoff"=unlist(ne[1:2]),
                                "Above cutoff"=unlist(po[1:2])))
    } else {
        ret <- withr::with_seed(42,
                                data.frame((hatM(cbind(Dm, Dp))[1:2])))
        rownames(ret) <- c("Pooled")
    }
    ret
}
