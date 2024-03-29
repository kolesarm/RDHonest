## kernel constants
kernels <- c("uniform", "triangular", "epanechnikov")
orders <- 0:2
kernC <- expand.grid(kernel=kernels, order=orders, boundary=c(TRUE, FALSE),
                     stringsAsFactors = FALSE)
kernC$mu0 <- 1
kernC$mu1 <- c(1/2, 1/3, 3/8, rep(0, nrow(kernC)-3))
kernC$mu2 <- c(1/3, 1/6, 1/5, -1/6, -1/10, -11/95, 0, 0, 0,
               1/3, 1/6, 1/5, 1/3, 1/6, 1/5, rep(0, 3))
kernC$mu3 <- c(1/4, 1/10, 1/8, -1/5, -1/10, -16/133, 1/20, 1/35, 15/448,
               rep(0, 9))
kernC$mu4 <- c(1/5, 1/15, 3/35, -1/5, -3/35, -141/1330, 3/35, 3/70, 13/252,
               1/5, 1/15, 3/35, 1/5, 1/15, 3/35, -3/35, -19/490, -1/21)
kernC$nu0 <- c(1,   4/3, 6/5, 4, 24/5, 56832/12635, 9, 72/7, 9895/1008,
               1/2, 2/3, 3/5, 1/2, 2/3, 3/5, 9/8, 456/343, 5/4)
kernC$nu1 <- c(1/2, 1/3, 3/8, 1, 3/5, 1770/2527, 3/2, 6/7, 4105/4032,
               rep(0, 9))
kernC$nu2 <- c(1/3, 2/15, 6/35, 8/15, 6/35, 2976/12635, 27/35, 8/35, 325/1008,
               1/6, 1/15, 3/35, 1/6, 1/15, 3/35, 9/56, 106/1715, 25/308)
kernC$nu3 <- c(1/4, 1/15, 3/32, 2/5, 3/35, 324/2527, 81/140, 4/35, 2825/16128,
               rep(0, 9))
kernC$nu4 <- c(1/5,  4/105, 2/35, 12/35, 2/35, 12352/138985, 17/35, 4/55,
               16795/144144,
               1/10, 2/105, 1/35, 1/10, 2/105, 1/35, 23/280, 256/18865, 85/4004)


kernC$pi0 <- c(1, 1, 1, 5/3, 3/2, 100483/64125,
               1+12*sqrt(6)/5^2, 1+2/sqrt(5), 2.0051585,
               rep(1, 6),
               1.323790008, 1.205511004, 1.244526871)
kernC$pi1 <- c(1/2, 1/3, 3/8, 16/27, 3/8, 702464/1603125,
               36*sqrt(6)/5^3, 24*sqrt(5)/5^3, 0.50792888,
               1/2, 1/3, 3/8, 1/2, 1/3, 3/8,
               0.4875, 0.3115591197, 0.3603316352)
kernC$pi2 <- c(1/3, 1/6, 1/5, 59/162, 3/16, 16520549/72140625,
               558*sqrt(6) / 5^5, 12*sqrt(5) / 5^3, 0.26617935,
               1/3, 1/6, 1/5, 1/3, 1/6, 1/5,
               0.2788548019, 0.1398694518, 0.1717750173)
kernC$pi3 <- c(1/4, 1/10, 1/8, 113/405, 1/8, 235792912/1514953125,
               1782*sqrt(6)/5^6+1/20, 218*sqrt(5)/4375+1/35,
               0.17770885,
               1/4, 1/10, 1/8, 1/4, 1/10, 1/8,
               0.1975, 0.08443054296, 0.1067100958)
kernC$pi4 <- c(1/5, 1/15, 3/35, 857/3645, 3/32, 1792724653/15149531250,
               0.268931, 0.102656, 0.1321148,
               1/5, 1/15, 3/35, 1/5, 1/15, 3/35,
               0.1574198065, 0.06014088703, 0.07706619387)

## Constant for pointwise MSE optimal bandwidth, page 67 in Fan and Gijbels
## (1996), pMSE=((p+1)!^2*nu / (2*(p+1)*mu_{p+1}^2))^{1/(2*p+3)}
kernC$pMSE <- NA
kernC$pMSE[kernC$order==0] <- (1/2 * kernC$nu0[kernC$order==0] /
                                   kernC$mu1[kernC$order==0]^2)^(1/3)
kernC$pMSE[kernC$order==1] <- (kernC$nu0[kernC$order==1] /
                                   kernC$mu2[kernC$order==1]^2)^(1/5)
kernC$pMSE[kernC$order==2] <- (6 * kernC$nu0[kernC$order==2] /
                                   kernC$mu3[kernC$order==2]^2)^(1/7)

## Save it
usethis::use_data(kernC, overwrite=TRUE, internal=TRUE)
