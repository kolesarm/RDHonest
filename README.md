[![R-CMD-check](https://github.com/kolesarm/RDHonest/workflows/R-CMD-check/badge.svg)](https://github.com/kolesarm/RDHonest/actions) [![Coverage status](https://codecov.io/gh/kolesarm/RDHonest/branch/master/graph/badge.svg)](https://codecov.io/github/kolesarm/RDHonest?branch=master)

# RDHonest

This R package implements honest and efficient confidence intervals in fuzzy and
sharp regression discontinuity designs using procedures from [Armstrong and
Kolesár (2020)](https://doi.org/10.3982/QE1199)
([preprint](https://arxiv.org/abs/1606.01200)), [Armstrong and Kolesár
(2018)](https://doi.org/10.3982/ECTA14434)
([preprint](https://arxiv.org/abs/1511.06028)), and [Kolesár and Rothe
(2018)](https://doi.org/10.1257/aer.20160945)
([preprint](https://arxiv.org/abs/1606.04086)).

See vignette [RDHonest](doc/RDHonest.pdf) for description of the package
(available through `vignette("RDHonest")` once package is installed), and the
package [manual](doc/manual.pdf) for documentation of the package functions.

This software package is based upon work supported by the National Science
Foundation under grant numbers SES-1628939 (Armstrong) and SES-1628878
(Kolesár).

## Installation

you can get the current development version from GitHub:

``` r
if (!requireNamespace("remotes")) {
  install.packages("remotes")
}
remotes::install_github("kolesarm/RDHonest")
```


TODO:
   - Format output
   - effective observations -> lindeberg weights.
   - allow for clustering
   - allow for covariates; For ROT: take the bandwidth that's optimal without
     covariates, for example, use it to compute the initial gamma, and then run
     the global quartic using the covariate adjusted outcome.

S3method(print,NPRResults)
S3method(print,RDBMEresults)
S3method(print,RDSmoothnessBound)
export(CVb)
export(EqKern)
export(FRDData)
export(FRDHonest)
export(KernMoment)
export(LPPData)
export(LPPHonest)
export(NPRPrelimVar.fit)
export(NPR_MROT.fit)
export(NPRreg.fit)
export(RDData)
export(RDHonest)
export(RDHonestBME)
export(RDSmoothnessBound)
export(RDTEfficiencyBound)
export(RDTOpt.fit)
export(plot_RDscatter)


TODO
- make class constructors internal
- for optimal Taylor inference, two bandwidths are reported

CHANGES
- simplify output of CVb function. No longer accept vector alpha as argument
- remove print method for NPRBW
- make NPRHonest.fit and NPROptbw.fit internal
- remove option to have different bandwidths on each side of cutoff.
- do not return lff
- fix bug with LPPOptBW whereby arguments to NPROptBW.fit not passed correctly
- remove functions FRDOptBW, RDOptBW and LPPOptBw
