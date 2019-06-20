#' Optimal bandwidth selection in RD
#'
#' Basic computing engine called by \code{\link{RDOptBW}} used to find
#' optimal bandwidth
#' @param d object of class \code{"RDData"}
#' @template RDoptBW
#' @template RDclass
#' @template Kern
#' @template bwequal
#' @template RDseInitial
#' @return a list with the following elements
#'     \describe{
#'     \item{\code{hp}}{bandwidth for observations above cutoff}
#'
#'     \item{\code{hm}}{bandwidth for observations below cutoff, equal to \code{hp}
#'     unless \code{bw.equal==FALSE}}
#'
#'     \item{\code{sigma2m}, \code{sigma2p}}{estimate of conditional variance
#'      above and below cutoff, from \code{d}}
#'    }
#' @references{
#' \cite{Imbens, Guido, and Kalyanaraman, Karthik,
#' "Optimal bandwidth choice for the regression discontinuity estimator." The
#' Review of Economic Studies 79 (3): 933-959.}
#' }
#' @examples
#' ## Lee data
#' d <- RDData(lee08, cutoff=0)
#' RDOptBW.fit(d, M=0.1, opt.criterion="MSE")[c("hp", "hm")]
#' @export
RDOptBW.fit <- function(d, M, kern="triangular", opt.criterion, bw.equal=TRUE,
                         alpha=0.05, beta=0.8, sclass="H", order=1,
                         se.initial="IKEHW", T0=0) {

    ## First check if sigma2 is supplied
    if (is.null(d$sigma2) & (is.null(d$sigma2p) | is.null(d$sigma2m)))
        d <- NPRPrelimVar.fit(d, se.initial=se.initial)

    ## Objective function for optimizing bandwidth
    obj <- function(hp, hm) {
        r <- RDHonest.fit(d, M, kern, c(p=abs(hp), m=abs(hm)),
                          alpha=alpha, se.method="supplied.var",
                          sclass=sclass, order=order)
        if (opt.criterion=="OCI") {
            2*r$maxbias+r$sd*(stats::qnorm(1-alpha)+stats::qnorm(beta))
        } else if (opt.criterion=="MSE") {
            r$maxbias^2+r$sd^2
        } else if (opt.criterion=="FLCI") {
            r$hl
        } else {
            stop(sprintf("optimality criterion %s not yet implemented",
                         opt.criterion))
        }
    }

    if (bw.equal==FALSE) {
        r <- stats::optim(c(max(d$Xp)/2, max(abs(d$Xm)/2)),
                   function(h) obj(h[1], h[2]))
        if (r$convergence!=0) warning(r$message)
        hp <- abs(r$par[1])
        hm <- abs(r$par[2])
    } else {
        obj1 <- function(h) obj(h, h)
        hmin <- max(unique(d$Xp)[order+1], sort(unique(abs(d$Xm)))[order+1])
        hmax <- max(abs(c(d$Xp, d$Xm)))
        ## Optimize piecewise constant function using modification of golden
        ## search. In fact, the criterion may not be unimodal, so proceed with
        ## caution. (For triangular kernel, it appears unimodal)
        if (kern=="uniform") {
            supp <- sort(unique(c(d$Xp, abs(d$Xm))))
            hp <- hm <- gss(obj1, supp[supp>=hmin])
        } else {
            hm <- hp <- abs(stats::optimize(obj1, interval=c(hmin, hmax),
                                        tol=.Machine$double.eps^0.75)$minimum)
        }

    }

    list(hp=hp, hm=hm, sigma2p=d$sigma2p, sigma2m=d$sigma2m)
}
