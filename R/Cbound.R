#' Lower bound on smoothness constant M in RD designs
#'
#' Estimate a lower bound on smoothness constant M and provide a lower
#' confidence interval.
#'
#' @param d object of class \code{"RDData"}
#' @param s Number of support points that curvature estimates should average
#'     over
#' @param sclass Smoothness class, either \code{"T"} for Taylor or \code{"H"}
#'     for Hölder class.
#' @param alpha determines confidence level \code{1-alpha}.
#' @param separate If \code{TRUE}, report estimates separately for data above
#'     and below cutoff. If \code{FALSE}, report pooled estimates
#' @param multiple If \code{TRUE}, use multiple curvature estimates. If
#'     \code{FALSE}, use a single estimate using only observations closest to
#'     the cutoff.
#' @return Returns a list with the following elements
#'
#' \describe{

#' \item{\code{mu+,mu-}}{Lower bound of CI for observations above and below
#'       cutoff}
#'
#' \item{\code{Z+,Z-}}{Point estimate used for lower bound}
#'
#' \item{\code{sd+,sd-}}{Standard deviations of point estimates}
#' }
#' @references{
#'
#' \cite{Armstrong, Timothy B., and Michal Kolesár. 2018.
#' "Optimal Inference in a Class of Regression Models." Econometrica 86 (2):
#' 655–83.}
#'
#' \cite{Kolesár, Michal, and Christoph Rothe. 2018. "Inference in Regression
#' Discontinuity Designs with a Discrete Running Variable." American Economic
#' Review 108 (8): 2277–2304.}
#' }
#' @export
RDSmoothnessBound <- function(d, s, separate=TRUE, multiple=TRUE, alpha=0.05,
                              sclass="T") {
    ## First estimate variance
    if (is.null(d$sigma2p) | is.null(d$sigma2m))
        d <- NPRPrelimVar.fit(d, se.initial="nn")

    ## Curvature estimate based on jth set of three points closest to zero
    Dk <- function(Y, X, xu, s2, j) {
        I1 <- X >= xu[3*j*s-3*s+1]  & X <= xu[3*j*s-2*s]
        I2 <- X >= xu[3*j*s-2*s+1]  & X <= xu[3*j*s-s]
        I3 <- X >= xu[3*j*s-  s+1]  & X <= xu[3*j*s]

        lam <- (mean(X[I3])-mean(X[I2])) / (mean(X[I3])-mean(X[I1]))
        den <- if  (sclass=="T") {
                   (1-lam)*mean(X[I3]^2) + lam*mean(X[I1]^2) + mean(X[I2]^2)
               } else {
                   (1-lam)*mean(X[I3]^2) + lam*mean(X[I1]^2) - mean(X[I2]^2)
               }
        Del <- 2*(lam*mean(Y[I1])+(1-lam)*mean(Y[I3])-mean(Y[I2])) / den
        VD <- 4*(lam^2*mean(s2[I1])/sum(I1) +
                 (1-lam)^2*mean(s2[I3])/sum(I3) + mean(s2[I2])/sum(I2)) / den^2
        c(Del, sqrt(VD), mean(Y[I1]), mean(Y[I2]), mean(Y[I3]), range(X[I1]),
          range(X[I2]), range(X[I3]))
    }

    xp <- unique(d$Xp)
    xm <- sort(unique(abs(d$Xm)))

    Dpj <- function(j) Dk(d$Yp, d$Xp, xp, d$sigma2p, j)
    Dmj <- function(j) Dk(d$Ym, abs(d$Xm), xm, d$sigma2m, j)

    if (multiple==TRUE) {
        ## Positive Deltas
        Sp <- floor(length(xp)/(3*s))
        Sm <- floor(length(xm)/(3*s))
        if (min(Sp, Sm) ==0) stop("Value of s is too big")
    } else {
        Sp <- Sm <- 1
    }
    Dp <- sapply(seq_len(Sp), Dpj)
    Dm <- sapply(seq_len(Sm), Dmj)

    ## Critical value
    cv <- function(M, Z, sd, alpha) {
        if (ncol(Z)==1) {
            return(CVb(M/sd, alpha=alpha)$cv)
        } else {
            S <- Z+M*outer(rep(1, nrow(Z)), 1/sd)
            maxS <- abs(S[cbind(seq_len(nrow(Z)), max.col(S))])
            return(unname(stats::quantile(maxS, 1-alpha)))
        }
    }

    hatM <- function(D) {
        ts <- abs(D[1, ]/D[2, ])
        maxt <- D[, which.max(ts)]
        set.seed(42)
        Z <- matrix(stats::rnorm(10000*ncol(D)), ncol=ncol(D))

        if (max(ts) < cv(0, Z, D[2, ], 1/2)) {
            hatM <- lower <- 0
        } else {
            hatM <- FindZero(function(m)
                max(ts)-cv(m, Z, D[2, ], 1/2), negative=FALSE)
            if (max(ts) < cv(0, Z, D[2, ], alpha)) {
                lower <- 0
            } else {
                lower <- FindZero(function(m)
                    max(ts)-cv(m, Z, D[2, ], alpha), negative=FALSE)
            }
        }
        list(hatM=hatM, lower=lower, Delta=maxt[1], sdDelta=maxt[2],
             y1=maxt[3], y2=maxt[4], y3=maxt[5],
             I1=maxt[6:7], I2=maxt[8:9], I3=maxt[10:11])
    }

    if (separate==TRUE) {
        po <- hatM(Dp)
        ne <- hatM(Dm)
    } else {
        po <- ne <- hatM(cbind(Dm, Dp))
    }

    ret <- list(po=po, ne=ne)
    class(ret) <- "RDSmoothnessBound"
    ret
}

#' @export
print.RDSmoothnessBound <- function(x, digits = getOption("digits"), ...) {
    fmt <- function(x) format(x, digits=digits, width=digits+1)

    pr <- function(r) {
        cat("Estimate: ", fmt(r$hatM), ", Lower CI: [", fmt(r$lower),
            ", Inf)\n", sep="")
        cat("\nDelta: ", fmt(r$Delta), ", sd=", fmt(r$sdDelta), sep="")
        cat("\nE_n[f(x_1)]: ", fmt(r$y1), ", I1=[", fmt(r$I1[1]), ", ",
            fmt(r$I1[2]), "]\n", "E_n[f(x_2)]: ", fmt(r$y2), ", I2=[",
            fmt(r$I2[1]), ", ", fmt(r$I2[2]), "]\n", "E_n[f(x_3)]: ", fmt(r$y3),
            ", I3=[", fmt(r$I3[1]), ", ", fmt(r$I3[2]), "]\n", sep="")
    }
    if (x$po$hatM==x$ne$hatM) {
        cat("\nSmoothness bound estimate:\n")
        pr(x$po)
    } else {
        cat("\nSmoothness bound estimate using observations above cutoff:\n")
        pr(x$po)
        cat("\nSmoothness bound estimate using observations below cutoff:\n")
        pr(x$ne)
    }

    invisible(x)
}
