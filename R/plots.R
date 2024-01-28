#' Scatterplot of binned raw observations
#'
#' Scatterplot of raw observations in which each point corresponds to an binned
#' average.
#'
#' @param cutoff specifies the RD cutoff for the running variable.
#' @template RDFormula
#' @param avg Number of observations to average over. If set to \code{Inf}, then
#'     take averages for each possible value of the running variable (convenient
#'     when the running variable is discrete).
#' @param xlab,ylab x- and y-axis labels
#' @param vert Draw a vertical line at cutoff?
#' @param propdotsize If \code{TRUE}, then size of points is proportional to
#'     number of observations that the point averages over (useful when
#'     \code{avg=Inf}). Otherwise the size of points is constant.
#' @return An object of class \code{"ggplot"}, a scatterplot the binned raw
#'     observations.
#' @examples
#' RDScatter(log(earnings)~yearat14, data=cghs, cutoff=1947,
#'                avg=Inf, propdotsize=TRUE)
#' @export
RDScatter <- function(formula, data, subset, cutoff=0, na.action, avg=10,
                      xlab=NULL, ylab=NULL, vert=TRUE, propdotsize=FALSE) {
    if (!requireNamespace("ggplot2", quietly = TRUE))
        stop("This function requires the ggplot2 package", call. = FALSE)
    ## construct model frame
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    d <- NPRData(mf, cutoff, "SRD")

    if (avg==Inf) {
        bd <- data.frame(x=unique(d$X),
                         y=unname(stats::coef(stats::lm(d$Y~0+as.factor(d$X)))),
                         count=as.vector(table(d$X)))
    } else {
        np <- sum(d$p)
        nm <- sum(d$m)
        ## don't recycle
        maxp <- (np %/% avg) * avg
        maxm <- (nm %/% avg) * avg
        bd <- data.frame(x=c(colMeans(matrix(d$X[d$m][1:maxm], nrow=avg)),
                             colMeans(matrix(d$X[d$p][1:maxp], nrow=avg))),
                         y=c(colMeans(matrix(d$Y[d$m][1:maxm], nrow=avg)),
                             colMeans(matrix(d$Y[d$p][1:maxp], nrow=avg))))
        ## if there is a remainder, add it
        if (maxm+1<=nm)
            bd <- rbind(bd, data.frame(x=mean(d$X[d$m][(maxm+1):nm]),
                                       y=mean(d$Y[d$m][(maxm+1):nm])))
        if (maxp+1<=np)
            bd <- rbind(bd, data.frame(x=mean(d$X[d$p][(maxp+1):np]),
                                       y=mean(d$Y[d$p][(maxp+1):np])))
        bd$count <- avg
    }
    bd$x <- bd$x+d$orig.cutoff
    if (propdotsize) {
        p <- ggplot2::ggplot()+
            ggplot2::geom_point(ggplot2::aes(x=bd$x, y=bd$y, size=bd$count))
    } else {
        p <- ggplot2::ggplot()+ggplot2::geom_point(ggplot2::aes(x=bd$x, y=bd$y))
    }

    p <- p + ggplot2::theme(legend.position = "none")
    if (!is.null(xlab)) p <- p + ggplot2::xlab(xlab)
    if (!is.null(ylab)) p <- p + ggplot2::ylab(ylab)
    if (vert) p <- p + ggplot2::geom_vline(xintercept=d$orig.cutoff,
                                           linetype="dotted")
    p
}
