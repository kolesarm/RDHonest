# RDHonest 1.0.1

## Minor improvements and fixes

- Use covariate-adjusted outcome to compute nearest-neighbor variance estimator
- Drop collinear covariates automatically instead of throwing an error


# RDHonest 1.0.0

## New Features

- The function `RDHonest` computes estimates and confidence intervals for the
  regression discontinuity (RD) parameter in sharp and fuzzy designs. It
  supports covariates, clustering, and weighting. Confidence intervals are
  honest (or bias-aware), with critical values computed using the `CVb`
  function. Worst-case bias of the estimator is computed under either the Taylor
  or Hölder smoothness class.
- `RDHonestBME` computes confidence intervals in sharp RD designs with discrete
  covariates under the assumption assumption that the conditional mean lies in
  the "bounded misspecification error" class of functions, as considered in
  [Kolesár and Rothe (2018)](https://doi.org/10.1257/aer.20160945).
- Support for plotting the data is provided by the function `RDScatter`
- The function `RDSmoothnessBound` computes a lower bound on the smoothness
  constant `M`, used as a parameter by `RDHonest` to calculate the worst-case
  bias of the estimator
- The function `RDTEfficiencyBound` calculates efficiency of minimax one-sided
  CIs at constant functions, or efficiency of two-sided fixed-length CIs at
  constant functions under second-order Taylor smoothness class.
