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
               orig.cutoff=cutoff, var.names=names(mf$model[1:2]))
    df$sigma2m <- d$"(sigma2)"[X<0]
    df$sigma2p <- d$"(sigma2)"[X>=0]

    structure(df, class="RDData")
}


#' check class of object
#' @keywords internal
CheckClass <- function(x, class)
    if(!inherits(x, class)) stop(paste0("Object ", deparse(substitute(x)),
                                        " needs to be class ", class, "!"))
