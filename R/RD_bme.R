#' Formula for local polynomial regression
#'
#' @param order order of local polynomial
#' @return Regression formula for local polynomial regression of order
#'     \code{order}
#' @keywords internal
RDlpformula <- function(order) {
        f1 <- if (order>0) {
                  f1 <- paste(sapply(1:order, function(p)
                      paste0("I(x^", p, ")")), collapse="+")
                  paste0("(", f1, ")*I(x>=0)")
              } else  {
                  paste0("I(x>=0)")
              }
        paste("y ~ ", f1)
}


#' BME method to compute honest CIs in RD designs with discrete regressors
#'
#' Computes honest CIs for local linear regression with uniform kernel as in
#' Kolesar and Rothe
#'
#' @template RDFormula
#' @template RDBW
#' @param alpha determines confidence level, \code{1-alpha}
#' @param order Order of local regression 1 for linear, 2 for quadratic.
#' @param regformula Explicitly specify regression formula as alternative to
#'     local linear regression, with \code{y} and \code{x} denoting the outcome
#'     and the running variable, and cutoff is normalized to \code{0}. Local
#'     linear regression (\code{order=1}) is equivalent to
#'     \code{regformula="y~x*I(x>0)"}. Inference is done the \code{order+2}th
#'     element of the design matrix
#' @importFrom stats model.matrix resid
#' @examples
#'
#' RDHonestBME(log(cghs$earnings)~yearat14, data=cghs, hp=3,
#'             order=1, cutoff=1947)
#' ## Equivalent to
#' RDHonestBME(log(cghs$earnings)~yearat14, data=cghs, hp=3,
#'             cutoff=1947, regformula="y~x*I(x>0)")
#' @export
RDHonestBME <- function(formula, data, subset, cutoff=0, na.action, hp=Inf,
                        hm=hp, alpha=0.05, order=0, regformula) {

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
    ind <- (x <= hp) & (x >= -hm)
    x <- x[ind]
    y <- stats::model.response(mf, "numeric")[ind]

    ## Count effective support points
    support <- sort(unique(x))
    G <- length(support)
    G.m <- length(support[support>=0])

    ## Estimate actual and dummied out model, and calculate delta
    m1 <- lm(regformula)
    m2 <- lm(y ~ 0+I(as.factor(x)))
    delta <- coef(m2)-stats::predict(m1, newdata=data.frame(x=support))

    ## Compute joint VCOV matrix of deltas and tau
    e2 <- rep(0, G+length(coef(m1)))
    e2[order + 2] <- 1                        # inference on (p+2)th element
    df <- data.frame(x=support, y=rep(0, length(support)))

    ## Compute Q^{-1} manually so that sandwich package is not needed
    Q1inv <- chol2inv(qr(m1)$qr[1L:m1$rank, 1L:m1$rank, drop = FALSE])
    Q2inv <- chol2inv(qr(m2)$qr[1L:m2$rank, 1L:m2$rank, drop = FALSE])
    v.m1m2 <- stats::var(cbind((model.matrix(m1)*resid(m1)) %*% Q1inv,
    (model.matrix(m2)*resid(m2)) %*% Q2inv)) * length(y)

    aa <- rbind(cbind(-model.matrix(regformula, data=df), diag(nrow=G)), e2)
    vdt <- aa %*% v.m1m2 %*% t(aa)

    ## All possible combinations of s_+, s- and g_+, g-
    gr <- as.matrix(expand.grid(1:G.m, (G.m+1):G, c(-1, 1), c(-1, 1)))
    selvec <- matrix(0, nrow=nrow(gr), ncol=ncol(vdt))
    selvec[cbind(1:nrow(selvec), gr[, 1])] <- gr[, 3]
    selvec[cbind(1:nrow(selvec), gr[, 2])] <- gr[, 4]
    selvec[, ncol(selvec)] <- 1

    se <- sqrt(rowSums((selvec %*% vdt) * selvec))
    dev <- drop(selvec[, -ncol(selvec)] %*% delta)

    ## Upper and lower CIs
    CI.l <- coef(m1)[order+2]+dev-qnorm(0.975)*se
    CI.u <- coef(m1)[order+2]+dev+qnorm(0.975)*se

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
    cat("(", fmt(r$CI[1]), ", ", fmt(r$CI[2]), ")\n", sep="")

    invisible(x)
}
