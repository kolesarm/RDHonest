#' Compute lower bound on C
#' @keywords internal
CboundXY <- function(Y, X, a, alpha, sigma2, sclass="T") {
    Ik <- function(A, k) a[k]<= A & A < a[k+1]
    E <- function(A, k) sum(A*Ik(X, k)) / sum(Ik(X, k))
    lam <- (E(X, 3)-E(X, 2)) / (E(X, 3)-E(X, 1))
    den <- if (sclass=="T") {
               lam*E(X^2, 1)+(1-lam)*E(X^2, 3)+E(X^2, 3)^2
           } else {
               stop("Don't know to compute the bound for this smoothness class")
           }
    Za <- (lam*E(Y, 1)+(1-lam)*E(Y, 3)-E(Y, 2)) / den
    sd <- lam^2*E(sigma2, 1)/sum(Ik(X, 1)) + (1-lam)^2*E(sigma2, 3) /
        sum(Ik(X, 3)) + E(sigma2, 2)/sum(Ik(X, 2))
    sd <- sqrt(sd) / den

    if (abs(Za/sd) < CVb(0, alpha=alpha)$cv) {
        hatmu <- 0
    } else {
        hatmu <- FindZero(function(m) abs(Za/sd)-CVb(m/sd, alpha=alpha)$cv,
                          negative=FALSE)
    }

    list(hatmu=hatmu, Zh=Za, sd=sd)
}


#' Lower bound on smoothness constant M
#'
#' Estimate a lower bound on smoothness constant M and provide a lower
#' confidence interval.
#'
#' @param d RDData object
#' @param ap,am Vectors specifying intervals to take averages over
#' @param sclass Smoothness class, either \code{"T"} for Taylor or
#'     \code{"H"} for Hölder class.
#' @param alpha determines confidence level \code{1-alpha}.
#' @param s with of a bin
#' @return Returns a list with the following elements
#'
#' \describe{

#' \item{\code{mu+,mu-}}{Lower bound of CI for observations above and below
#'       cutoff}}
#'
#' \item{\code{Z+,Z-}}{Point estimate used for lower bound}}
#'
#' \item{\code{sd+,sd-}}{Standard deviations of point estimates}}
#' }
#' @export
RDMbound <- function(d, s, ap, am, alpha=0.5, sclass="T") {
    if (missing(ap) | missing(am)) {
        ap <- c(0, d$Xp[s+1], d$Xp[2*s+1], d$Xp[3*s+1])
        nm <- length(d$Xm)
        am <- c(0, abs(d$Xm)[nm-s], abs(d$Xm)[nm-2*s], abs(d$Xm)[nm-3*s])
    }

    rp <- CboundXY(d$Yp, d$Xp, ap, alpha, mean(d$sigma2p), sclass)
    rm <-  CboundXY(d$Ym, abs(d$Xm), am, alpha, mean(d$sigma2m), sclass)
    names(rp) <- c("mu+", "Z+", "sd+")
    names(rm) <- c("mu-", "Z-", "sd-")
    r <- unlist(c(rp, rm))
    ## Convert C to M
    2*r
}
