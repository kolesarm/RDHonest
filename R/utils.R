## Class Constructor
NPRData <- function(d, cutoff, class, fo) {
    ## Drop intercepts
    rhs <- function(j) stats::model.matrix(stats::formula(fo, lhs=0, rhs=j), d)

    X <- rhs(1)[, -1, drop=FALSE]-cutoff
    msg <- paste0("Formula syntax:\n'outcome ~ running_variable |",
                  " covariates' for sharp RD,\n'outcome | treatment ~ ",
                  "running_variable | covariates' for fuzzy RD",
                  ",\n'outcome ~ running_variable' for inference at a point.")
    if (ncol(X) != 1)
        stop(paste0("Single running variable required.\n", msg))
    if (is.unsorted(X)) {
        idx <- sort(X, index.return=TRUE)$ix
        X <- X[idx, , drop=FALSE]
        d <- d[idx, ]
    }
    Y <- Formula::model.part(fo, data = d, lhs = 1)
    if (length(fo)[1] > 1)
        Y <- cbind(Y, Formula::model.part(fo, data = d, lhs = 2))

    df <- structure(list(X=X, Y=as.matrix(Y), orig.cutoff=cutoff,
                         p=drop(X>=0), m=drop(X<0),
                         clusterid=d$"(clusterid)", w=stats::model.weights(d)),
                    class=class)
    if (ncol(df$Y)!=length(fo)[1])
        stop(paste0("Single outcome and treatment variable required.", msg))
    if (is.null(df$w)) df$w <- rep(1L, NROW(X))
    if (length(fo)[2] > 1) df$covs <- rhs(2)[, -1, drop=FALSE]

    if (length(fo)[1] == 1) {
        df$sigma2 <- d$"(sigmaY2)"
    } else if (!is.null(d$"(sigmaY2)")) {
        df$sigma2 <- cbind(d$"(sigmaY2)", d$"(sigmaYD)", d$"(sigmaYD)",
                           d$"(sigmaD2)")
    }
    df
}


## Find interval containing zero of a function, then find the zero
## Search an interval for a root of \code{f},
## @param f function whose root we're looking for
## @param ival upper endpoint of initial interval in which to search
## @param negative logical: should the lower endpoint be \code{1/ival} (if the
##     root is guaranteed to be positive), or \code{-ival}?
FindZero <- function(f, ival=1.1, negative=TRUE) {
    minval <- function(ival) if (negative==TRUE) -ival else min(1/ival, 1e-3)

    while (sign(f(ival))==sign(f(minval(ival))))
        ival <- 2*ival
    stats::uniroot(f, c(minval(ival), ival), tol=.Machine$double.eps^0.75)$root
}


## Modified golden section for unimodal piecewise constant function
gss <- function(f, xs) {
    gr <- (sqrt(5) + 1) / 2
    a <- 1
    b <- length(xs)
    c <- round(b - (b - a) / gr)
    d <- round(a + (b - a) / gr)

    while (b - a > 100) {
        if (f(xs[c]) < f(xs[d])) {
            b <- d
        } else {
            a <- c
        }
        # recompute both c and d
        c <- round(b - (b - a) / gr)
        d <- round(a + (b - a) / gr)
    }

    supp <- xs[a:b]
    supp[which.min(vapply(supp, f, numeric(1)))]
}
