## ---- include=FALSE, cache=FALSE-----------------------------------------
library("knitr")
knitr::opts_knit$set(self.contained = FALSE)
knitr::opts_chunk$set(tidy = TRUE, collapse=TRUE, comment = "#>",
                      tidy.opts=list(blank=FALSE, width.cutoff=55))

## ---- fig.width=4.5, fig.height=3.5, fig.cap="Lee (2008) data"-----------
library("RDHonest")
## transform data to an RDdata object
dt <- RDData(lee08, cutoff = 0)
## plot 25-bin averages in a window equal to 50 around the cutoff, see Figure 1
plot_RDscatter(dt, avg=25, window = 50) +
   xlab("Margin of victory") + ylab("Vote share in next election")

## ------------------------------------------------------------------------
## Usual critical value
CVb(0, alpha=0.05)
## Tabulate critical values for different significance levels
## when bias-sd ratio equals 1/4
knitr::kable(CVb(1/4, alpha=c(0.01, 0.05, 0.1)))

## ------------------------------------------------------------------------
RDHonest(voteshare ~ margin, data=lee08, kern="uniform", M=0.1, hp=10, sclass="T")
RDHonest(voteshare ~ margin, data=lee08, kern="uniform", M=0.1, hp=10, sclass="H")

## ------------------------------------------------------------------------
RDHonest(voteshare ~ margin, data=lee08, kern="uniform", M=0.1, opt.criterion="MSE", sclass="T")
RDHonest(voteshare ~ margin, data=lee08, kern="uniform", M=0.1, opt.criterion="MSE", sclass="H")

## ------------------------------------------------------------------------
RDOptBW(voteshare ~ margin, data=lee08, kern="uniform", M=0.1, opt.criterion="MSE", sclass="T")
RDOptBW(voteshare ~ margin, data=lee08, kern="uniform", M=0.1, opt.criterion="MSE", sclass="H")

