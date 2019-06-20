#' Honest inference in RD
#'
#' Calculate estimators and one- and two-sided CIs based on local polynomial
#' estimator in RD under second-order Taylor or Hölder smoothness class. If
#' \code{kern="optimal"}, calculate optimal estimators under second-order Taylor
#' smoothness class.
#'
#' The bandwidth is calculated to be optimal for a given performance criterion,
#' as specified by \code{opt.criterion}. For local polynomial estimators, this
#' optimal bandwidth is calculated using the function \code{\link{RDOptBW}}.
#' Alternatively, for local polynomial estimators, the bandwidths above and
#' below the cutoff can be specified by \code{h}.
#'
#' @template FRDFormula
#' @template RDse
#' @template RDoptBW
#' @template RDBW
#' @template RDclass
#' @template Kern
#' @template bwequal
#' @template RDseInitial
#' @return Returns an object of class \code{"RDResults"}. The function
#'     \code{print} can be used to obtain and print a summary of the results. An
#'     object of class \code{"RDResults"} is a list containing the following
#'     components
#'
#'     \describe{
#'   \item{\code{estimate}}{Point estimate. This estimate is MSE-optimal if
#'                   \code{opt.criterion="MSE"}}
#'
#'   \item{lff}{Least favorable function, only relevant for optimal estimator
#'              under Taylor class.}
#'
#'   \item{\code{maxbias}}{Maximum bias of \code{estimate}}
#'
#'   \item{\code{sd}}{Standard deviation of estimate}
#'
#'   \item{\code{lower}, \code{upper}}{Lower (upper) end-point of a one-sided CI
#'         based on \code{estimate}. This CI is optimal if
#'         \code{opt.criterion=="OCI"}}
#'
#'   \item{\code{hl}}{Half-length of a two-sided CI based on \code{estimate}, so
#'             that the CI is given by \code{c(estimate-hl, estimate+hl)}. The
#'             CI is optimal if \code{opt.criterion="FLCI"}}
#'
#'   \item{\code{eff.obs}}{Effective number of observations used by
#'             \code{estimate}}
#'
#'   \item{\code{hp}, \code{hm}}{Bandwidths used}
#'
#'   \item{\code{naive}}{Coverage of CI that ignores bias and uses
#'                \code{qnorm(1-alpha/2)} as critical value}
#'
#'   \item{\code{call}}{the matched call}
#'
#' }
#' @seealso \code{\link{RDOptBW}}
#' @references{
#' \cite{Imbens, Guido, and Kalyanaraman, Karthik,
#' "Optimal bandwidth choice for the regression discontinuity estimator." The
#' Review of Economic Studies 79 (3): 933-959.}
#' }
#' @examples
#'
#' # Lee dataset
#' RDHonest(voteshare ~ margin, data = lee08, kern = "uniform", M = 0.1,
#'          h = 10, sclass = "T")
#' @export
FRDHonest <- function(formula, data, subset, cutoff=0, M, kern="triangular",
                     na.action, opt.criterion, bw.equal=TRUE, hp, hm=hp,
                     se.method="nn", alpha=0.05, beta=0.8, J=3, sclass="H",
                     order=1, se.initial="IKEHW") {

    ## construct model frame
    cl <- mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    if (!requireNamespace("Formula", quietly = TRUE)) {
        stop("Please install the Formula package;
              it's needed for this function to work",
             call. = FALSE)
    }


    formula <- Formula::as.Formula(formula)
    stopifnot(length(formula)[1] == 1L, length(formula)[2] == 2)
    mf$formula <- formula

    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    d <- FRDData(mf, cutoff)

    if (!missing(hp)) {
        ret <- FRDHonest.fit(d, M, kern, hp, hm, alpha=alpha,
                            se.method=se.method, J=J, sclass=sclass,
                            order=order, se.initial=se.initial)
    } else {
        ret <- FRDHonest.fit(d, M, kern, opt.criterion=opt.criterion,
                             bw.equal=bw.equal, alpha=alpha, beta=beta,
                             se.method=se.method, J=J,
                             sclass=sclass, order=order, se.initial=se.initial)
    }

    ret$call <- cl
    ret$na.action <- attr(mf, "na.action")

    ret
}


#' Optimal Bandwidth Selection in Regression Discontinuity
#'
#' Estimate bandwidth for sharp RD based on local polynomial regression that
#' optimizes either maximum mean squared error, or length or quantiles of excess
#' length of a honest CI under second order Hölder or Taylor class.
#'
#' @template RDFormula
#' @template RDoptBW
#' @template RDclass
#' @template Kern
#' @template bwequal
#' @template RDseInitial
#' @return Returns an object of class \code{"RDBW"}. The function \code{print}
#'     can be used to obtain and print a summary of the results. An object of
#'     class \code{"RDBW"} is a list containing the following components:
#'
#'     \describe{
#'     \item{\code{hp}}{bandwidth for observations above cutoff}
#'
#'     \item{\code{hm}}{bandwidth for observations below cutoff, equal to \code{hp}
#'     unless \code{bw.equal==FALSE}}
#'
#'     \item{\code{sigma2m}, \code{sigma2p}}{estimate of conditional variance
#'     just above and just below cutoff, \eqn{\sigma^2_+(0)} and
#'     \eqn{\sigma^2_{-}(0)}}
#'
#'    \item{\code{f0}}{estimate of density of running variable at cutoff, if
#'    bandwidth computed using asymptotic method}
#'
#'    \item{\code{call}}{the matched call}
#'
#'    \item{\code{na.action}}{(where relevant) information on handling of missing
#'    data.}
#'
#'    }
#' @seealso \code{\link{RDHonest}}
#' @references{
#' \cite{Imbens, Guido, and Kalyanaraman, Karthik,
#' "Optimal bandwidth choice for the regression discontinuity estimator." The
#' Review of Economic Studies 79 (3): 933-959.}
#' }
#' @examples
#'
#' ## Use Lee dataset
#' RDOptBW(voteshare ~ margin, data = lee08, kern = "uniform",
#'         M = 0.1, opt.criterion = "MSE", sclass = "H")
#' @export
FRDOptBW <- function(formula, data, subset, cutoff=0, M, kern="triangular",
                    na.action, opt.criterion, bw.equal=TRUE,
                    alpha=0.05, beta=0.8, sclass="H", order=1,
                    se.initial="IKEHW") {

    ## construct model frame
    cl <- mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    d <- RDData(mf, cutoff)

    ret <- RDOptBW.fit(d, M, kern, opt.criterion, bw.equal, alpha, beta, sclass,
                       order, se.initial=se.initial)
    class(ret) <- "RDBW"
    ret$call <- cl
    ret$na.action <- attr(mf, "na.action")

    ret
}

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
FRDOptBW.fit <- function(d, M, kern="triangular", opt.criterion,
                        bw.equal=TRUE, alpha=0.05, beta=0.8,
                        sclass="H", order=1, se.initial="IKEHW", T0=0) {

    ## First check if sigma2 is supplied
    if (is.null(d$sigma2p) | is.null(d$sigma2m))
        d <- FRDPrelimVar(d, se.initial=se.initial)

    ## Objective function for optimizing bandwidth
    obj <- function(hp, hm) {
        r <- FRDHonest.fit(d, M, kern, c(p=abs(hp), m=abs(hm)),
                           alpha=alpha, se.method="supplied.var",
                           sclass=sclass, order=order, T0=T0)
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

#' @export
print.FRDResults <- function(x, digits = getOption("digits"), ...) {
    if (!is.null(x$call))
        cat("Call:\n", deparse(x$call), "\n\n", sep = "", fill=TRUE)
    bw <- if (class(x$lff) !="RDLFFunction")
              "Bandwidth" else "Smoothing parameter"

    cat("Inference by se.method:\n")
    r <- as.data.frame(x[c("estimate", "maxbias", "sd",
                           "lower", "upper", "hl")])
    names(r)[1:3] <- c("Estimate", "Maximum Bias", "Std. Error")
    r$name <- rownames(r)
    print.data.frame(as.data.frame(r[, 1:3]), digits=digits)
    cat("\nConfidence intervals:\n")
    fmt <- function(x) format(x, digits=digits, width=digits+1)

    for (j in seq_len(nrow(r)))
        cat(format(r[j, 7], width=5), " (", fmt(r[j, 1]-r[j, 6]),
            ", ", fmt(r[j, 1]+r[j, 6]), "), (", fmt(r[j, 4]),
            ", Inf), (-Inf, ", fmt(r[j, 5]), ")\n", sep="")

    cat("\n", bw, " below cutoff: ", format(x$hm, digits=digits), sep="")
    cat("\n", bw, " above cutoff: ", format(x$hp, digits=digits), sep="")
    if (x$hm==x$hp) {
        cat(" (", bw, "s are the same)\n", sep="")
    } else {
        cat(" (", bw, "s are different)\n", sep="")
    }
    cat("Number of effective observations:",
        format(x$eff.obs, digits=digits), "\n")

    invisible(x)
}
