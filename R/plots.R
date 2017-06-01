#' Scatterplot of binned raw observations
#'
#' Scatterplot of raw observations in which each point corresponds to an binned
#' average.
#' @param d RDdata object
#' @param avg is number of observations to average over. If set to \code{Inf},
#'     then take averages for each possible value of the running variable
#'     (convenient when the running variable is discrete)
#' @param xlab,ylab x- and y-axis labels
#' @param window Width of a window around cutoff to which the graph should be
#'     restricted. If not specified, full data range will be plotted
#' @param vert Draw a vertical line at cutoff?
#' @param propdotsize If \code{TRUE}, then size of points is proportional to
#'     numer of observations that the point averages over (useful when
#'     \code{avg=Inf}). Otherwise the size of points is constant.
#' @export
plot_RDscatter <- function(d, avg=10, xlab=NULL, ylab=NULL,
                           window=NULL, vert=TRUE, propdotsize=FALSE) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Please install ggplot2 package;
              it's needed for this function to work",
             call. = FALSE)
    }

    ## RDData is sorted
    if (!is.null(window)) {
        d$Yp <- d$Yp[d$Xp<=window]
        d$Ym <- d$Ym[d$Xm>=-window]
        d$Xp <- d$Xp[d$Xp<=window]
        d$Xm <- d$Xm[d$Xm>=-window]
    }

    if (avg==Inf) {
        x <- c(d$Xm, d$Xp)
        y <- c(d$Ym, d$Yp)
        bd <- data.frame(x=unique(x), y=unname(coef(lm(y~0+as.factor(x)))),
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
        p <- ggplot2::qplot(x=x, y=y, data=bd, size=count)
    } else {
        p <- ggplot2::qplot(x=x, y=y, data=bd)
    }

    p <- p + ggplot2::theme_classic() +
        ggplot2::theme(legend.position = "none", axis.line.x=ggplot2::
                       element_line(color="black", size=0.2, linetype="solid"),
                       axis.line.y=ggplot2::
                       element_line(color="black", size=0.2, linetype="solid"))
    if(!is.null(xlab)) p <- p + ggplot2::xlab(xlab)
    if(!is.null(ylab)) p <- p + ggplot2::ylab(ylab)
    if(vert) p <- p + ggplot2::geom_vline(xintercept=d$orig.cutoff,
                                          linetype="dotted")

    p
}
