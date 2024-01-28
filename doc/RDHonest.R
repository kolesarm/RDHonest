## ----include=FALSE, cache=FALSE-----------------------------------------------
library("knitr")
knitr::opts_knit$set(self.contained = FALSE)
knitr::opts_chunk$set(tidy = TRUE, collapse=TRUE, comment = "#>",
                      tidy.opts=list(blank=FALSE, width.cutoff=55))

## -----------------------------------------------------------------------------
library("RDHonest")

## ----fig.width=4.5, fig.height=3.5, fig.cap="Lee (2008) data"-----------------
## plot 25-bin averages in for observations 50 at most points
## away from the cutoff. See Figure 1.
RDScatter(voteshare~margin, data=lee08, subset=abs(lee08$margin)<=50,
          avg=50, propdotsize=FALSE, xlab="Margin of victory",
          ylab="Vote share in next election")

## ----fig.width=4.5, fig.height=3.5, fig.cap="Oreopoulos (2006) data"----------
## see Figure 2
f2 <- RDScatter(I(log(earnings))~yearat14, data=cghs, cutoff=1947,
                avg=Inf, xlab="Year aged 14", ylab="Log earnings",
                propdotsize=TRUE)
## Adjust size of dots if they are too big
f2 + ggplot2::scale_size_area(max_size = 4)

## -----------------------------------------------------------------------------
CVb(0, alpha=0.05) ## Usual critical value
CVb(1/2, alpha=0.05)

## Tabulate critical values for different bias levels
CVb(0:5, alpha=0.1)

## -----------------------------------------------------------------------------
RDHonest(voteshare~margin, data=lee08, kern="uniform", M=0.1, h=10, sclass="T")
RDHonest(voteshare~margin, data=lee08, kern="uniform", M=0.1, h=10, sclass="H")

## -----------------------------------------------------------------------------
RDHonest(voteshare ~ margin, data=lee08, kern="triangular",
         M=0.1, opt.criterion="MSE", sclass="H")
## Choose bws optimal for length of CI
RDHonest(voteshare ~ margin, data=lee08, kern="triangular", M=0.1,
         opt.criterion="FLCI", sclass="H")

## -----------------------------------------------------------------------------
## Replicate Table 2, column (10)
RDHonest(log(earnings) ~ yearat14, cutoff=1947,
         data=cghs, kern="uniform", M=0.04, opt.criterion="FLCI", sclass="H")
## Triangular kernel generally gives tigher CIs
RDHonest(log(earnings) ~ yearat14, cutoff=1947,
         data=cghs, kern="triangular", M=0.04, opt.criterion="FLCI", sclass="H")

## -----------------------------------------------------------------------------
## Replicate Table 2, column (6), run local linear regression (order=1)
## with a uniform kernel (other kernels are not yet implemented)
RDHonestBME(log(earnings) ~ yearat14, cutoff=1947,
            data=cghs, h=3, order=1)

## -----------------------------------------------------------------------------
## Data-driven choice of M
RDHonest(voteshare ~ margin, data=lee08, kern="uniform", sclass="H",
         opt.criterion="MSE")

## -----------------------------------------------------------------------------
r1 <- RDHonest(voteshare ~ margin, data=lee08, kern="optimal", M=0.1,
               opt.criterion="FLCI",
               se.method="nn")$coefficients
r2 <- RDHonest(voteshare ~ margin, data=lee08, kern="triangular", M=0.1,
               opt.criterion="FLCI", se.method="nn",
               sclass="T")$coefficients
r1$conf.high-r1$conf.low
r2$conf.high-r2$conf.low

## -----------------------------------------------------------------------------
## Add variance estimate to the Lee (2008) data so that the RDSmoothnessBound
## function doesn't have to compute them each time
## dl <- PrelimVar(dl, se.initial="nn")

### Only use three point-average for averages of a 100 points closest to cutoff,
### and report results separately for points above and below cutoff
## RDSmoothnessBound(dl, s=100, separate=TRUE, multiple=FALSE, sclass="T")

### Pool estimates based on observations below and above cutoff, and
### use three-point averages over the entire support of the running variable
## RDSmoothnessBound(dl, s=100, separate=FALSE, multiple=TRUE, sclass="H")

## -----------------------------------------------------------------------------
d <- cghs
## Make 20 groups based on observation number
d$mod <- seq_along(d$yearat14) %% 20
## Make cells defined as intersection of group and year
d$cell <- d$mod/100+d$yearat14
## Data with cell averages
dd <- data.frame()
for (j in unique(d$cell)){
    dd <- rbind(dd, data.frame(y=mean(log(d$earnings)[d$cell==j]),
                               x=mean(d$yearat14[d$cell==j]),
                               weights=length(d$yearat14[d$cell==j])))
}

## -----------------------------------------------------------------------------
RDHonest(log(earnings)~yearat14, cutoff=1947, h=5, data=cghs, M=0.1,
         se.method="nn")
RDHonest(y~x, cutoff=1947, weights=weights, h=5, data=dd, M=0.1,
         se.method="nn")

## -----------------------------------------------------------------------------
## Assumes first column in the data frame corresponds to outcome,
##  second to the treatment variable, and third to the running variable
## Outcome here is log of non-durables consumption
## dr <- FRDData(cbind(logf=log(rcp[, 6]), rcp[, c(3, 2)]), cutoff=0)

## -----------------------------------------------------------------------------
## Initial estimate of treatment effect for optimal bandwidth calculations
r <- RDHonest(log(cn) ~ retired | elig_year, data=rcp, kern="triangular",
              M=c(0.001, 0.002), opt.criterion="MSE", sclass="H", T0=0)
## Use it to compute optimal bandwidth
RDHonest(log(cn) ~ retired | elig_year, data=rcp, kern="triangular",
         M=c(0.001, 0.002), opt.criterion="MSE", sclass="H",
         T0=r$coefficients$estimate)

## -----------------------------------------------------------------------------
## Data-driven choice of M
RDHonest(log(cn) ~ retired | elig_year, data=rcp, kern="triangular",
         opt.criterion="MSE", sclass="H", T0=r$coefficients$estimate)

## -----------------------------------------------------------------------------
## Transform data, specify we're interested in inference at x0=20,
## and drop observations below cutoff
RDHonest(voteshare ~ margin, data=lee08, subset=margin>0,
         cutoff=20, kern="uniform",
         opt.criterion="MSE", sclass="H", point.inference=TRUE)

