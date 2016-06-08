#' Scatterplot of binned raw observations
#'
#' Scatterplot of raw observations in which each point corresponds to an binned
#' average.
#' @param d RDdata object
#' @param avg is number of observations to average over
#' @param xlab,ylab x- and y-axis labels
#' @param window Width of a window around cutoff to which the graph should be
#'     restricted. If not specified, full data range will be plotted
#' @param vert Draw a vertical line at cutoff?
#' @export
plot_RDscatter <- function(d, avg=10, xlab=NULL, ylab=NULL,
                           window=NULL, vert=TRUE) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Please install ggplot2 package;
              it's needed for this function to work",
             call. = FALSE)
    }

    ## Sort data
    if (is.unsorted(d$Xp) | is.unsorted(d$Xm)) {
        s <- sort(d$Xp, index.return=TRUE)
        d$Yp <- d$Yp[s$ix]
        d$Xp <- s$x
        s <- sort(d$Xm, index.return=TRUE)
        d$Ym <- d$Ym[s$ix]
        d$Xm <- s$x
    }
    if (!is.null(window)) {
        d$Yp <- d$Yp[d$Xp<=window]
        d$Ym <- d$Ym[d$Xm>=-window]
        d$Xp <- d$Xp[d$Xp<=window]
        d$Xm <- d$Xm[d$Xm>=-window]
    }

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
    ## if there is a remainded, add it
    if (maxm+1<=nm)
        bd <- rbind(bd, data.frame(x=mean(d$Xm[(maxm+1):nm]),
                                   y=mean(d$Ym[(maxm+1):nm])))
    if (maxp+1<=np)
        bd <- rbind(bd, data.frame(x=mean(d$Xp[(maxp+1):np]),
                                   y=mean(d$Yp[(maxp+1):np])))

    p <- ggplot2::qplot(x=x, y=y, data=bd) + ggplot2::theme_classic() +
        ggplot2::theme(axis.line.x=ggplot2::
                       element_line(color="black", size=0.2, linetype="solid"),
                       axis.line.y=ggplot2::
                       element_line(color="black", size=0.2, linetype="solid"))
    if(!is.null(xlab)) p <- p + ggplot2::xlab(xlab)
    if(!is.null(ylab)) p <- p + ggplot2::ylab(ylab)
    if(vert) p <- p + ggplot2::geom_vline(xintercept=0, linetype="dotted")

    p
}


## PlotBWRange <- function(Y, X, bws=NULL, kernel="triangular", basebw=NULL){
##     ## plot bandwidth over a given range

##     df <- data.frame(bw=bws, theta=NA, se=NA, n.plus=NA, n.minus=NA)

##     for (row in 1:nrow(df)) {
##         r <- RDreg(Y, X, bw=df$bw[row], kernel=kernel, degree=1, which.se="nn")
##         df[row, -1] <- c(r$theta, r$se["1nn"], sum(X<=df$bw[row] & X>=0),
##                          sum(X>=-df$bw[row] & X<0))
##     }

##     equiv.kernel <- paste("l", kernel, sep="")
##     df$cv <- LookUpCV(h.ratio=max(bws)/min(bws), kernel=equiv.kernel)

##     p <- ggplot(data=df) +  geom_line(aes(x=bw, y=theta)) + theme_bw() +
##         geom_line(aes(x=bw, y=theta+1.96*se), linetype="dotted") +
##         geom_line(aes(x=bw, y=theta-1.96*se), linetype="dotted") +
##         geom_line(aes(x=bw, y=(theta+cv*se)), linetype="dashed", color="blue") +
##         geom_line(aes(x=bw, y=(theta-cv*se)), linetype="dashed", color="blue")

##     if(!is.null(basebw))
##         p <- p + geom_vline(xintercept=basebw, linetype="longdash")

##     return(list(plot=p, table=df))
## }
