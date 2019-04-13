[![Travis build status](https://travis-ci.org/kolesarm/RDHonest.svg?branch=master)](https://travis-ci.org/kolesarm/RDHonest) [![Coverage status](https://codecov.io/gh/kolesarm/RDHonest/branch/master/graph/badge.svg)](https://codecov.io/github/kolesarm/RDHonest?branch=master)

# RDHonest

Honest and efficient confidence intervals in regression discontinuity designs
using procedures from [Armstrong and Kolesár
(2018)](https://arxiv.org/abs/1606.01200), [Armstrong and Kolesár
(2018)](https://doi.org/10.3982/ECTA14434)
([preprint](https://arxiv.org/abs/1511.06028)), and [Kolesár and Rothe
(2018)](https://doi.org/10.1257/aer.20160945)
([preprint](https://arxiv.org/abs/1606.04086))

See vignette [RDHonest](doc/RDHonest.pdf) for description of the package
(available through `vignette("RDHonest")` once package is installed), and the
package [manual](doc/manual.pdf) for documentation of the package functions.

This software package is based upon work supported by the National Science
Foundation under grant numbers SES-1628939 (Armstrong) and SES-1628878
(Kolesár).

## Installation

You can install the package manually by downloading the source code here, or
using the function `install_github()` from the `devtools` package:

```
install.packages("devtools") ## if devtools package not installed
devtools::install_github("kolesarm/RDHonest")
```
