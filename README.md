[![R-CMD-check](https://github.com/kolesarm/RDHonest/workflows/R-CMD-check/badge.svg)](https://github.com/kolesarm/RDHonest/actions) [![Coverage status](https://codecov.io/gh/kolesarm/RDHonest/branch/master/graph/badge.svg)](https://app.codecov.io/github/kolesarm/RDHonest?branch=master) [![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/RDHonest)](https://cran.r-project.org/package=RDHonest) [![Download statistics](https://cranlogs.r-pkg.org/badges/RDHonest)](https://cran.r-project.org/package=RDHonest)

# RDHonest

This R package implements honest and efficient confidence intervals in fuzzy and
sharp regression discontinuity designs using procedures from [Armstrong and
Kolesár (2020)](https://doi.org/10.3982/QE1199)
([preprint](https://arxiv.org/abs/1606.01200)), [Armstrong and Kolesár
(2018)](https://doi.org/10.3982/ECTA14434)
([preprint](https://arxiv.org/abs/1511.06028)), and [Kolesár and Rothe
(2018)](https://doi.org/10.1257/aer.20160945)
([preprint](https://arxiv.org/abs/1606.04086)). See
[RDHonest-vStata](https://github.com/tbarmstr/RDHonest-vStata) for a Stata
version of this package.

See vignette [RDHonest](doc/RDHonest.pdf) for description of the package
(available through `vignette("RDHonest")` once package is installed), and the
package [manual](doc/manual.pdf) for documentation of the package functions.

This software package is based upon work supported by the National Science
Foundation under grant numbers SES-1628939 (Armstrong) and SES-1628878
(Kolesár).

## Installation

You can install the released version of `RDHonest` from
[CRAN](https://CRAN.R-project.org/package=RDHonest) with:

``` r
install.packages("RDHonest")
```


Alternatively, you can get the current development version from GitHub:

``` r
if (!requireNamespace("remotes")) {
  install.packages("remotes")
}
remotes::install_github("kolesarm/RDHonest")
```
