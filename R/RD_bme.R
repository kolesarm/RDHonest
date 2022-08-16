## Formula for local polynomial regression
RDlpformula <- function(order) {
    f1 <- if (order>0) {
              f <- vapply(seq_len(order),
                           function(p) paste0("I(x^", p, ")"),
                           character(1))
              paste0("(", paste(f, collapse="+"), ") * I(x>=0)")
          } else  {
              paste0("I(x>=0)")
          }
    paste0("y ~ ", f1)
}


#' Honest CIs in sharp RD with discrete regressors under BME function class
#'
#' Computes honest CIs for local polynomial regression with uniform kernel under
#' the assumption that the conditional mean lies in the bounded misspecification
#' error (BME) class of functions, as considered in Kolesár and Rothe (2018).
#' This class formalizes the notion that the fit of the chosen model is no worse
#' at the cutoff than elsewhere in the estimation window.
#'
#' @template RDFormula
#' @param cutoff specifies the RD cutoff in the running variable.
#' @param h bandwidth, a scalar parameter.
#' @param alpha determines confidence level, \eqn{1-\alpha}{1-alpha}
#' @param order Order of local regression \code{1} for linear, \code{2} for
#'     quadratic, etc.
#' @param regformula Explicitly specify regression formula to use instead of
#'     running a local polynomial regression, with \code{y} and \code{x}
#'     denoting the outcome and the running variable, and cutoff is normalized
#'     to \code{0}. Local linear regression (\code{order = 1}) is equivalent to
#'     \code{regformula = "y~x*I(x>0)"}. Inference is done on the
#'     \code{order+2}th element of the design matrix
#' @return An object of class \code{"RDResults"}. This is a list with at least
#'     the following elements:
#'
#'    \describe{
#'
#'    \item{\code{"coefficients"}}{Data frame containing estimation results,
#'    including point estimate, one- and two-sided confidence intervals, a bound
#'    on worst-case bias, bandwidth used, and the number of effective
#'    observations.}
#'
#'    \item{\code{"call"}}{The matched call.}
#'
#'    \item{\code{"na.action"}}{(If relevant) information on the special
#'    handling of \code{NA}s.}
#'
#' }
#' @examples
#' RDHonestBME(log(earnings)~yearat14, data=cghs, h=3,
#'             order=1, cutoff=1947)
#' ## Equivalent to
#' RDHonestBME(log(earnings)~yearat14, data=cghs, h=3,
#'             cutoff=1947, order=1, regformula="y~x*I(x>=0)")
#' @references{
#'
#' \cite{Michal Kolesár and Christoph Rothe. Inference in regression
#'       discontinuity designs with a discrete running variable. American
#'       Economic Review, 108(8):2277—-2304, August 2018.
#'       \doi{10.1257/aer.20160945.}}
#'
#' }
#' @export
RDHonestBME <- function(formula, data, subset, cutoff=0, na.action,
                        h=Inf, alpha=0.05, order=0, regformula) {
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
    ind <- (x <= h) & (x >= -h)
    x <- x[ind]
    y <- stats::model.response(mf, "numeric")[ind]

    ## Count effective support points
    support <- sort(unique(x))
    G <- length(support)
    G.m <- length(support[support<0])

    ## Estimate actual and dummied out model, and calculate delta
    m1 <- stats::lm(regformula)
    m2 <- stats::lm(y ~ 0+I(as.factor(x)))
    delta <- stats::coef(m2)-stats::predict(m1, newdata=data.frame(x=support))

    ## Compute joint VCOV matrix of deltas and tau
    ## Compute Q^{-1} manually so that sandwich package is not needed
    Q1inv <- chol2inv(qr(m1)$qr[1L:m1$rank, 1L:m1$rank, drop = FALSE])
    Q2inv <- chol2inv(qr(m2)$qr[1L:m2$rank, 1L:m2$rank, drop = FALSE])
    v.m1m2 <- length(y) * stats::var(cbind(
    (stats::model.matrix(m1)*stats::resid(m1)) %*% Q1inv,
    (stats::model.matrix(m2)*stats::resid(m2)) %*% Q2inv))

    df <- data.frame(x=support, y=rep(0, length(support)))
    e2 <- rep(0, G+length(stats::coef(m1)))
    e2[order + 2] <- 1                  # inference on (p+2)th element
    aa <- rbind(cbind(-stats::model.matrix(regformula, data=df), diag(nrow=G)),
                e2)
    vdt <- aa %*% v.m1m2 %*% t(aa)      # V(W) in paper, except order of m1 and
                                        # m2 swapped

    ## All possible combinations of g_-, g_+, s_-, s_+
    gr <- as.matrix(expand.grid(1:G.m, (G.m+1):G, c(-1, 1), c(-1, 1)))
    selvec <- matrix(0, nrow=nrow(gr), ncol=ncol(vdt))
    selvec[cbind(seq_len(nrow(selvec)), gr[, 1])] <- gr[, 3]
    selvec[cbind(seq_len(nrow(selvec)), gr[, 2])] <- gr[, 4]
    selvec[, ncol(selvec)] <- 1

    se <- sqrt(rowSums((selvec %*% vdt) * selvec))
    dev <- drop(selvec[, -ncol(selvec)] %*% delta)

    ## Upper and lower CIs
    CI.l <- stats::coef(m1)[order+2]+dev-stats::qnorm(1-alpha/2)*se
    CI.u <- stats::coef(m1)[order+2]+dev+stats::qnorm(1-alpha/2)*se
    ## Onesided
    OCI.l <- stats::coef(m1)[order+2]+dev-stats::qnorm(1-alpha)*se
    OCI.u <- stats::coef(m1)[order+2]+dev+stats::qnorm(1-alpha)*se

    l <- which.min(CI.l)
    u <- which.max(CI.u)
    coef <- data.frame(term="Sharp RD parameter",
                       estimate=unname(m1$coefficients[order +2]),
                       std.error=sqrt(vdt["e2", "e2"]),
                       maximum.bias=max(abs(c(dev[u], dev[l]))),
                       conf.low=CI.l[l], conf.high=CI.u[u],
                       conf.low.onesided=min(OCI.l),
                       conf.high.onesided=max(OCI.u), bandwidth=h,
                       eff.obs=length(x), cv=NA, alpha=alpha, method="BME",
                       kernel="uniform")
    ret <- list(coefficients=coef, call=cl, na.action=attr(mf, "na.action"))
    class(ret) <- "RDResults"

    ret
}
