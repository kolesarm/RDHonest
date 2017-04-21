tol <- .Machine$double.eps^0.75

#' Find interval containing zero of a function, then find the zero
#'
#' Same function as in NPRHonest
#'
#' Given function \code{f} find \code{x0} such that \code{f(x0)==0}
#' @param f function whose root we're looking for
#' @param ival upper endpoint of initial interval in which to search
#' @param negative logical: should the lower endpoint be \code{1/ival} (if the
#'     root is guaranteed to be positive), or \code{-ival}?
#' @keywords internal
FindZero <- function(f, ival=1.1, negative=TRUE) {
    minval <- function(ival) if (negative==TRUE) -ival else min(1/ival, 1e-3)

    while(sign(f(ival))==sign(f(minval(ival))))
            ival <- 2*ival
    stats::uniroot(f, c(minval(ival), ival), tol=tol)$root
}


#' Class Constructor for RDData
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
#'     \item{Yp}{Outcome vector for observationsabove cutoff}
#'
#'     \item{Xm}{Running variable for observations below cutoff}
#'
#'     \item{Xp}{Running variable for observationsabove cutoff}
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
#' ## TODO
#' @export
RDData <- function(d, cutoff) {

    X <- d[[2]] - cutoff
    df <- list(Ym=d[[1]][X<0], Yp=d[[1]][X>=0], Xm=X[X<0], Xp=X[X>=0],
               orig.cutoff=cutoff, var.names=names(d)[1:2])
    df$sigma2m <- d$"(sigma2)"[X<0]
    df$sigma2p <- d$"(sigma2)"[X>=0]

    structure(df, class="RDData")
}


#' check class of object
#' @keywords internal
CheckClass <- function(x, class)
    if(!inherits(x, class)) stop(paste0("Object ", deparse(substitute(x)),
                                        " needs to be class ", class, "!"))


#' Split function into k bits and optimize on each bit in case not convex
#' @keywords internal
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
