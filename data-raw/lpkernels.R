## ---- echo=TRUE, results=FALSE------------------------------------------------
library("RDHonest")
EqKern("uniform", boundary = TRUE, order = 2)(0.5)
# Equivalent call
EqKern(function(u) u <= 1, boundary = TRUE, order = 2)(0.5)

## -----------------------------------------------------------------------------
## mu_1, should be 0
KernMoment(function(u) 4-6*u, moment = 1, boundary = TRUE, type = "raw")
## nu_1, should be 0
KernMoment(function(u) 4-6*u, moment = 1, boundary = TRUE, type = "raw2")
## pi_1, should be 16/27
KernMoment(function(u) 4-6*u, moment = 1, boundary = TRUE, type = "absolute")

