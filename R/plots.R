#' Scatterplot of binned raw observations
#'
#' Scatterplot of raw observations in which each point corresponds to an binned
#' average.
#'
#' @param formula object of class \code{"formula"} (or one that can be coerced
#'     to that  class) of the form \code{y ~ x}
#' @param data optional data frame, list or environment (or object coercible by
#'     \code{as.data.frame} to a data frame) containing the outcome and running
#'     variables in the model. If not found in \code{data}, the variables are
#'     taken from \code{environment(formula)}, typically the environment from
#'     which the function is called.
#' @param subset optional vector specifying a subset of observations to be used
#' @param cutoff specifies the RD cutoff in the running variable.
#' @param avg Number of observations to average over. If set to \code{Inf}, then
#'     take averages for each possible value of the running variable (convenient
#'     when the running variable is discrete).
#' @param xlab,ylab x- and y-axis labels
#' @param vert Draw a vertical line at cutoff?
#' @param propdotsize If \code{TRUE}, then size of points is proportional to
#'     number of observations that the point averages over (useful when
#'     \code{avg=Inf}). Otherwise the size of points is constant.
#' @return A \code{"ggplot"} object.
#' @examples
#' plot_RDscatter(I(log(earnings))~yearat14, data=cghs, cutoff=1947,
#'                avg=Inf, propdotsize=TRUE)
#' @export
plot_RDscatter <- function(formula, data, subset, cutoff=0, avg=10,
                           xlab=NULL, ylab=NULL, vert=TRUE, propdotsize=FALSE) {
    if (!requireNamespace("ggplot2", quietly = TRUE))
        stop("This function requires the ggplot2 package to be installed",
             call. = FALSE)
    ## construct model frame
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    d <- RDData(mf, cutoff)

    if (avg==Inf) {
        x <- c(d$Xm, d$Xp)
        y <- c(d$Ym, d$Yp)
        bd <- data.frame(x=unique(x),
                         y=unname(stats::coef(stats::lm(y~0+as.factor(x)))),
                         count=as.vector(table(x)))
    } else {
        np <- length(d$Yp)
        nm <- length(d$Ym)
        ## don't recycle
        maxp <- (np %/% avg) * avg
        maxm <- (nm %/% avg) * avg
        bd <- data.frame(
            x=c(colMeans(matrix(d$Xm[1:maxm], nrow=avg)),
                colMeans(matrix(d$Xp[1:maxp], nrow=avg))),
            y=c(colMeans(matrix(d$Ym[1:maxm], nrow=avg)),
                colMeans(matrix(d$Yp[1:maxp], nrow=avg))))
        ## if there is a remainder, add it
        if (maxm+1<=nm)
            bd <- rbind(bd, data.frame(x=mean(d$Xm[(maxm+1):nm]),
                                       y=mean(d$Ym[(maxm+1):nm])))
        if (maxp+1<=np)
            bd <- rbind(bd, data.frame(x=mean(d$Xp[(maxp+1):np]),
                                       y=mean(d$Yp[(maxp+1):np])))
        bd$count <- avg
    }
    bd$x <- bd$x+d$orig.cutoff

    if (propdotsize) {
        p <- ggplot2::qplot(x=bd$x, y=bd$y, size=bd$count)
    } else {
        p <- ggplot2::qplot(x=x, y=y, data=bd)
    }

    p <- p + ggplot2::theme(legend.position = "none")
    if(!is.null(xlab)) p <- p + ggplot2::xlab(xlab)
    if(!is.null(ylab)) p <- p + ggplot2::ylab(ylab)
    if(vert) p <- p + ggplot2::geom_vline(xintercept=d$orig.cutoff,
                                          linetype="dotted")

    p
}
