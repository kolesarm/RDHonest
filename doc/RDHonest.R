## ----include=FALSE, cache=FALSE-----------------------------------------------
library("knitr")
knitr::opts_knit$set(self.contained = FALSE)
knitr::opts_chunk$set(tidy = TRUE, collapse=TRUE, comment = "#>",
                      tidy.opts=list(blank=FALSE, width.cutoff=55))

## ----fig.width=4.5, fig.height=3.5, fig.cap="Lee (2008) data"-----------------
library("RDHonest")
## plot 50-bin averages in for observations 50 at most points
## away from the cutoff. See Figure 1.
RDScatter(voteshare~margin, data=lee08, subset=abs(lee08$margin)<=50,
          avg=50, propdotsize=FALSE, xlab="Margin of victory",
          ylab="Vote share in next election")

## ----fig.width=4.5, fig.height=3.5, fig.cap="Oreopoulos (2006) data"----------
## see Figure 2
f2 <- RDScatter(log(earnings)~yearat14, data=cghs, cutoff=1947,
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
r0 <- RDHonest(voteshare~margin, data=lee08, kern="triangular", M=0.1, h=8)
print(r0)

## -----------------------------------------------------------------------------
RDHonest(voteshare ~ margin, data=lee08, kern="triangular",
         M=0.1, opt.criterion="MSE")
## Choose bws optimal for length of CI
RDHonest(voteshare ~ margin, data=lee08, kern="triangular", M=0.1,
         opt.criterion="FLCI")

## -----------------------------------------------------------------------------
## Data-driven choice of M
RDHonest(voteshare ~ margin, data=lee08)

## -----------------------------------------------------------------------------
## Replicate Table 2, column (10) in Kolesar and Rothe (2018)
RDHonest(log(earnings) ~ yearat14, cutoff=1947,
         data=cghs, kern="uniform", M=0.04, opt.criterion="FLCI", sclass="H")

## -----------------------------------------------------------------------------
## Replicate Table 2, column (6), run local linear regression (order=1)
## with a uniform kernel (other kernels are not implemented for RDHonestBME)
RDHonestBME(log(earnings) ~ yearat14, cutoff=1947,
            data=cghs, h=3, order=1)

## -----------------------------------------------------------------------------
## Initial estimate of treatment effect for optimal bandwidth calculations
r <- RDHonest(log(cn) | retired ~ elig_year, data=rcp, kern="triangular",
              M=c(0.001, 0.002), opt.criterion="MSE", sclass="H", T0=0)
## Use it to compute optimal bandwidth
RDHonest(log(cn) | retired ~ elig_year, data=rcp, kern="triangular",
         M=c(0.001, 0.002), opt.criterion="MSE", sclass="H",
         T0=r$coefficients$estimate)

## -----------------------------------------------------------------------------
## Data-driven choice of M
RDHonest(log(cn) | retired ~ elig_year, data=rcp, kern="triangular",
         opt.criterion="MSE", sclass="H", T0=r$coefficients$estimate)

## -----------------------------------------------------------------------------
## No covariates
rn <- RDHonest(mortHS ~ povrate, data=headst)
## Use Percent attending school aged 14-17, urban, black,
## and their interaction as covariates.
rc <- RDHonest(mortHS ~ povrate | urban*black + sch1417, data=headst)
rn
rc

## -----------------------------------------------------------------------------
ci_len <- c(rc$coefficients$conf.high-rc$coefficients$conf.low,
            rn$coefficients$conf.high-rn$coefficients$conf.low)
100 * (1 - ci_len[1]/ci_len[2])

## -----------------------------------------------------------------------------
RDHonest(log(cn) | retired ~ elig_year | education, data=rcp,
         T0=r$coefficients$estimate)

## ----fig.width=4.5, fig.height=3.5, fig.cap="Battistin et al (2009) data"-----
## see Figure 3
f3 <- RDScatter(log(cn)~elig_year, data=rcp, cutoff=0, avg=Inf,
                xlab="Years to eligibility",
                ylab="Log consumption of non-durables", propdotsize=TRUE,
                subset=abs(elig_year)<10)
## Adjust size of dots if they are too big
f3 + ggplot2::scale_size_area(max_size = 5)

## -----------------------------------------------------------------------------
dd <- data.frame()
## Collapse data by running variable
for (j in unique(cghs$yearat14)) {
    ix <- cghs$yearat14==j
    df <- data.frame(y=mean(log(cghs$earnings[ix])), x=j,
                     weights=sum(ix),
                     sigma2=var(log(cghs$earnings[ix]))/sum(ix))
    dd <- rbind(dd, df)
}

## -----------------------------------------------------------------------------
s0 <- RDHonest(log(earnings)~yearat14, cutoff=1947, data=cghs)
## keep same bandwidth
s1 <- RDHonest(y~x, cutoff=1947, data=dd, weights=weights,
               sigmaY2=sigma2, se.method="supplied.var",
               h=s0$coefficients$bandwidth)
## Results are identical:
s0
s1

## -----------------------------------------------------------------------------
r0 <- RDHonest(log(cn)|retired~elig_year, data=rcp, h=7)
dd <- data.frame(x=sort(unique(rcp$elig_year)), y=NA, d=NA, weights=NA,
                 sig11=NA, sig12=NA, sig21=NA, sig22=NA)
for (j in seq_len(NROW(dd))) {
    ix <- rcp$elig_year==dd$x[j]
    Y <- cbind(log(rcp$cn[ix]), rcp$retired[ix])
    dd[j, -1] <- c(colMeans(Y), sum(ix), as.vector(var(Y))/sum(ix))
}
r1 <- RDHonest(y|d~x, data=dd, weights=weights,
               sigmaY2=sig11, T0=0, sigmaYD=sig21,
               sigmaD2=sig22, h=7,
               se.method="supplied.var")
## Outputs match up to numerical precision
max(abs(r0$coefficients[2:11]-r1$coefficients[2:11]))

## -----------------------------------------------------------------------------
## make fake clusters
set.seed(42)
clusterid <- sample(1:50, NROW(lee08), replace=TRUE)
sc <- RDHonest(voteshare~margin, data=lee08, se.method="EHW",
               clusterid=clusterid, M=0.14, h=7)
## Since clusters are unrelated to outcomes, not clustering
## should yield similar standard errors
sn <- RDHonest(voteshare~margin, data=lee08, se.method="EHW",
               M=0.14, h=7)
sc
sn

## -----------------------------------------------------------------------------
r1 <- RDHonest(voteshare ~ margin, data=lee08, M=0.1, se.method="nn")
### Only use three point-average for averages of a 100 points closest to cutoff,
### and report results separately for points above and below cutoff
RDSmoothnessBound(r1, s=100, separate=TRUE, multiple=FALSE, sclass="T")

## Pool estimates based on observations below and above cutoff, and
## use three-point averages over the entire support of the running variable
RDSmoothnessBound(r1, s=100, separate=FALSE, multiple=TRUE, sclass="H")

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
## Specify we're interested in inference at x0=20,
## and drop observations below cutoff
RDHonest(voteshare ~ margin, data=lee08, subset=margin>0,
         cutoff=20, kern="uniform",
         opt.criterion="MSE", sclass="H", point.inference=TRUE)

