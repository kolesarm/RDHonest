tol <- .Machine$double.eps^0.75

## Find interval containing zero of a function, then find the zero

## Given function \code{f} find \code{x0} such that \code{f(x0)==0}
## @param f function whose root we're looking for
## @param ival upper endpoint of initial interval in which to search
## @param negative logical: should the lower endpoint be \code{1/ival} (if the
##     root is guaranteed to be positive), or \code{-ival}?
FindZero <- function(f, ival=1.1, negative=TRUE) {
    minval <- function(ival) if (negative==TRUE) -ival else min(1/ival, 1e-3)

    while(sign(f(ival))==sign(f(minval(ival))))
            ival <- 2*ival
    stats::uniroot(f, c(minval(ival), ival), tol=tol)$root
}


#' Class Constructor for \code{"RDData"}
#'
#' Convert data to standardized format for use with low-level functions. If the
#' cutoff for treatment is non-zero, shift the running variable so that cutoff
#' is at zero.
#'
#' @param d data frame with first column corresponding to outcome variable,
#'     second column corresponding to running variable and optionally a column
#'     called \code{"(sigma2)"} that corresponds to the conditional variance of
#'     the outcome (or an estimate of the conditional variance)
#' @param cutoff specifies the RD cutoff in the running variable.
#' @return An object of class \code{"RDData"}, which is a list containing the
#'     following components:
#'
#'     \describe{
#'
#'     \item{Ym}{Outcome vector for observations below cutoff}
#'
#'     \item{Yp}{Outcome vector for observations above cutoff}
#'
#'     \item{Xm}{Running variable for observations below cutoff}
#'
#'     \item{Xp}{Running variable for observations above cutoff}
#'
#'     \item{sigma2m}{Conditional variance of outcome for observations below
#'     cutoff}
#'
#'     \item{sigma2p}{Conditional variance of outcome for observations above
#'     cutoff}
#'
#'     \item{orig.cutoff}{Original cutoff}
#'
#'     item{var.names}{Names of outcome and running variable in supplied data frame}
#'
#'     }
#' @examples
#'
#' ## Transform Lee data
#' d <- RDData(lee08, cutoff=0)
#' @export
RDData <- function(d, cutoff) {

    X <- d[[2]] - cutoff
    df <- list(Ym=d[[1]][X<0], Yp=d[[1]][X>=0], Xm=X[X<0], Xp=X[X>=0],
               orig.cutoff=cutoff, var.names=names(d)[1:2])
    df$sigma2m <- d$"(sigma2)"[X<0]
    df$sigma2p <- d$"(sigma2)"[X>=0]

    # Sort data
    s <- sort(df$Xp, index.return=TRUE)
    df$Yp <- df$Yp[s$ix]
    df$Xp <- s$x
    s <- sort(df$Xm, index.return=TRUE)
    df$Ym <- df$Ym[s$ix]
    df$Xm <- s$x

    structure(df, class="RDData")
}

#' Class Constructor for \code{"LPPData"}
#'
#' Convert data to standardized format for use with low-level functions. If the
#' point for which to do inference is non-zero, shift the independent variable so
#' that it is at zero.
#'
#' @param d data frame with first column corresponding to outcome variable,
#'     second column corresponding to independent variable and optionally a column
#'     called \code{"(sigma2)"} that corresponds to the conditional variance of
#'     the outcome (or an estimate of the conditional variance)
#' @param point specifies the point at which to calculate conditional mean
#' @return An object of class \code{"LPPData"}, which is a list containing the
#'     following components:
#'
#'     \describe{
#'
#'     \item{Y}{Outcome vector}
#'
#'     \item{X}{Independent variable}
#'
#'     \item{sigma2}{Conditional variance of outcome}
#'
#'     \item{orig.point}{Original point}
#'
#'     item{var.names}{Names of outcome and independent variable in supplied data frame}
#'
#'     }
#' @examples
#'
#' ## Transform Lee data
#' d <- RDData(lee08[lee08$margin>0, ], cutoff=0)
#' @export
LPPData <- function(d, point) {

    X <- d[[2]] - point
    df <- list(Y=d[[1]], X=X,
               orig.point=point, var.names=names(d)[1:2])
    df$sigma2 <- d$"(sigma2)"

    # Sort data
    s <- sort(df$X, index.return=TRUE)
    df$Y <- df$Y[s$ix]
    df$X <- s$x

    structure(df, class="LPPData")
}


## check class of object
CheckClass <- function(x, class)
    if(!inherits(x, class)) stop(paste0("Object ", deparse(substitute(x)),
                                        " needs to be class ", class, "!"))


## Split function into k bits and optimize on each bit in case not convex
CarefulOptim <- function(f, interval, k=10) {
    ## intervals
    s <- seq(interval[1], interval[2], length.out=k+1)
    s <- matrix(c(s[-length(s)], s[-1]), ncol=2)

    obj <- rep(0, k)
    arg <- rep(0, k)
    for (j in 1:k){
        r <- stats::optimize(f, s[j, ])
        arg[j] <- r$minimum
        obj[j] <- r$objective
    }
    jopt <- which.min(obj)
    list(objective=obj[jopt], minimum=arg[jopt])
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
