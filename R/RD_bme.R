## Formula for local polynomial regression
## @param order order of local polynomial
## @return Regression formula for local polynomial regression of order
##     \code{order}
RDlpformula <- function(order) {
        f1 <- if (order>0) {
                  f1 <- paste(vapply(seq_len(order), function(p)
                      paste0("I(x^", p, ")"), character(1)), collapse="+")
                  paste0("(", f1, ")*I(x>=0)")
              } else  {
                  paste0("I(x>=0)")
              }
        paste("y ~ ", f1)
}


#' CIs in sharp RD with discrete regressors under bounded misspecification error
#' class
#'
#' Computes honest CIs for local linear regression with uniform kernel under the
#' bounded misspecification error class of functions, as considered in Kolesár
#' and Rothe (2018)
#'
#' The parameter \code{weights} is ignored, it is only included to keep a
#' unified interface with \code{\link{RDHonest}}.
#'
#' @template RDFormula
#' @template RDBW
#' @param alpha determines confidence level, \eqn{1-\alpha}{1-alpha}
#' @param order Order of local regression \code{1} for linear, \code{2} for
#'     quadratic.
#' @param regformula Explicitly specify regression formula to use instead of
#'     running a local linear regression, with \code{y} and \code{x} denoting
#'     the outcome and the running variable, and cutoff is normalized to
#'     \code{0}. Local linear regression (\code{order = 1}) is equivalent to
#'     \code{regformula = "y~x*I(x>0)"}. Inference is done on the
#'     \code{order+2}th element of the design matrix
#' @examples
#' RDHonestBME(log(cghs$earnings)~yearat14, data=cghs, h=3,
#'             order=1, cutoff=1947)
#' ## Equivalent to
#' RDHonestBME(log(cghs$earnings)~yearat14, data=cghs, h=3,
#'             cutoff=1947, order=1, regformula="y~x*I(x>=0)")
#' @references{
#'
#' \cite{Kolesár, Michal, and Christoph Rothe. 2018. "Inference in Regression
#' Discontinuity Designs with a Discrete Running Variable." American Economic
#' Review 108 (8): 2277–2304.}
#'
#' }
#' @export
RDHonestBME <- function(formula, data, subset, weights, cutoff=0, na.action,
                        h=Inf, alpha=0.05, order=0, regformula) {
    if (length(h)==1) h <- c(m=h, p=h)
    ## construct model frame
    cl <- mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    if(missing(regformula))
        regformula <- RDlpformula(order)
    regformula <- stats::as.formula(regformula)

    x <- mf[, 2]-cutoff
    ## drop observations outside bw
    ind <- (x <= h["p"]) & (x >= -h["m"])
    x <- x[ind]
    y <- stats::model.response(mf, "numeric")[ind]

    ## Count effective support points
    support <- sort(unique(x))
    G <- length(support)
    G.m <- length(support[support>=0])

    ## Estimate actual and dummied out model, and calculate delta
    m1 <- stats::lm(regformula)
    m2 <- stats::lm(y ~ 0+I(as.factor(x)))
    delta <- stats::coef(m2)-stats::predict(m1, newdata=data.frame(x=support))

    ## Compute joint VCOV matrix of deltas and tau
    ## Compute Q^{-1} manually so that sandwich package is not needed
    Q1inv <- chol2inv(qr(m1)$qr[1L:m1$rank, 1L:m1$rank, drop = FALSE])
    Q2inv <- chol2inv(qr(m2)$qr[1L:m2$rank, 1L:m2$rank, drop = FALSE])
    v.m1m2 <- stats::var(cbind((
        stats::model.matrix(m1)*stats::resid(m1)) %*% Q1inv,
    (stats::model.matrix(m2)*stats::resid(m2)) %*% Q2inv)) * length(y)


    df <- data.frame(x=support, y=rep(0, length(support)))
    e2 <- rep(0, G+length(stats::coef(m1)))
    e2[order + 2] <- 1                  # inference on (p+2)th element
    aa <- rbind(cbind(-stats::model.matrix(
                                  regformula, data=df), diag(nrow=G)), e2)
    vdt <- aa %*% v.m1m2 %*% t(aa)      # V(W) in paper, except order of m1 and
                                        # m2 swapped

    ## All possible combinations of s_+, s- and g_+, g-
    gr <- as.matrix(expand.grid(1:G.m, (G.m+1):G, c(-1, 1), c(-1, 1)))
    selvec <- matrix(0, nrow=nrow(gr), ncol=ncol(vdt))
    selvec[cbind(seq_len(nrow(selvec)), gr[, 1])] <- gr[, 3]
    selvec[cbind(seq_len(nrow(selvec)), gr[, 2])] <- gr[, 4]
    selvec[, ncol(selvec)] <- 1

    se <- sqrt(rowSums((selvec %*% vdt) * selvec))
    dev <- drop(selvec[, -ncol(selvec)] %*% delta)

    ## Upper and lower CIs
    CI.l <- stats::coef(m1)[order+2]+dev-stats::qnorm(0.975)*se
    CI.u <- stats::coef(m1)[order+2]+dev+stats::qnorm(0.975)*se

    l <- gr[which.min(CI.l), ]
    deltamin <- c(support[l[1:2]], dev[which.min(CI.l)],
                  se[which.min(CI.l)])
    u <- gr[which.max(CI.u), ]
    deltamax <- c(support[u[1:2]], dev[which.max(CI.u)],
                  se[which.max(CI.u)])
    names(deltamin) <- names(deltamax) <- c("x_{-}", "x_{+}", "delta", "se")

    ret <- list(CI=unname(c(min(CI.l), max(CI.u))),
                delta.l=deltamin, delta.u=deltamax, call=cl,
                na.action=attr(mf, "na.action"), regformula=regformula,
                x=x, y=y)
    class(ret) <- "RDBMEresults"
    ret
}

#' @export
print.RDBMEresults <- function(x, digits = getOption("digits"), ...) {
    if (!is.null(x$call))
        cat("Call:\n", deparse(x$call), "\n\n", sep = "", fill=TRUE)
    fmt <- function(x) format(x, digits=digits, width=digits+1)

    cat("Confidence intervals:\n")
    cat("(", fmt(x$CI[1]), ", ", fmt(x$CI[2]), ")\n", sep="")

    invisible(x)
}
