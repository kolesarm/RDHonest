context("Test RD")

test_that("Test class constructor sorting", {
    d0 <- NPRData(rcp[, c(6, 3, 2)], cutoff=3, "FRD")
    d1 <- NPRData(rcp[sort(rcp$elig_year,
                           index.return=TRUE)$ix, c(6, 3, 2)], cutoff=3, "FRD")
    expect_identical(d1, d0)

    d0 <- NPRData(rcp[, c(6, 2)], cutoff=3, "SRD")
    d1 <- NPRData(rcp[sort(rcp$elig_year,
                           index.return=TRUE)$ix, c(6, 2)], cutoff=3, "SRD")
    expect_identical(d1, d0)

    d0 <- NPRData(rcp[, c(6, 2)], cutoff=-3, "IP")
    d1 <- NPRData(rcp[sort(rcp$elig_year,
                           index.return=TRUE)$ix, c(6, 2)], cutoff=-3, "IP")
    expect_identical(d1, d0)

    ## Now test that I() works
    expect_equal(RDHonest(voteshare ~ margin, data=lee08,
                          M=0, h=2)$coefficients$estimate,
                 -RDHonest(voteshare ~ I(-margin), data=lee08,
                           M=0, h=2)$coefficients$estimate)
    expect_equal(RDHonest(cn~retired | elig_year, data=rcp, cutoff=0,
                          M=c(1, 0.1), h=3)$coefficients$estimate,
                 RDHonest(cn~retired | I(2*elig_year), data=rcp, cutoff=0,
                          M=c(1, 0.1), h=6)$coefficients$estimate)
})

test_that("IK bandwidth calculations", {
    ## Test IK bandwidth in Lee data, IK Table 1
    d <- NPRData(lee08, cutoff=0, "SRD")

    dig <- getOption("digits")
    options(digits=8)
    r1 <- capture.output(IKBW(d, verbose=TRUE))
    expect_equal(r1[c(2, 4, 5, 6, 8)],
                 c(" h1:  14.44507 ", " f(0):  0.0089622411 ",
                   " sigma^2_{+}(0):  12.024424 ^2",
                   " sigma^2_{+}(0): 10.472064 ^2 ",
                   " h_{2, +}: 60.513312 h_{2, -}: 60.993358 "))
    options(digits=dig)

    expect_equal(IKBW(d), 29.3872649956)

    r <- NPReg(d, IKBW(d, kern="uniform"), "uniform")
    expect_equal(r$estimate, 8.0770003749)
    d <- PrelimVar(NPRData(lee08, cutoff=0, "SRD"), se.initial="EHW")
    expect_equal(sqrt(mean(d$sigma2[d$p])), 12.58183131)
    expect_equal(sqrt(mean(d$sigma2[d$m])), 10.79067278)
})

test_that("Plots", {
    if (requireNamespace("ggplot2", quietly = TRUE)) {
        expect_silent(invisible(RDScatter(earnings~yearat14, data=cghs,
                                          cutoff=1947, avg=Inf,
                                          propdotsize=TRUE)))
        expect_silent(invisible(RDScatter(voteshare~margin, data=lee08,
                                          subset=abs(lee08$margin)<=50, avg=50,
                                          propdotsize=FALSE, xlab="margin",
                                          ylab="effect")))
    }
})


test_that("Honest inference in Lee and LM data",  {

    ## Replicate 1606.01200v2, except we no longer provide SilvermanNN
    r1 <- RDHonest(mortHS ~ povrate, data=headst, kern="uniform", h=18, M=0.1,
                   se.method="EHW")
    r2 <- RDHonest(mortHS ~ povrate, data=headst, kern="uniform", h=18, M=0.1,
                   se.method="nn")
    ## Should match regression
    rl <- lm(mortHS ~ povrate*I(povrate>=0), data=headst,
             subset=abs(povrate)<=18)
    XX <- model.matrix(rl)
    meat <- crossprod(XX, rl$residuals^2*XX)
    vl <- (solve(crossprod(XX)) %*% meat %*% solve(crossprod(XX)))[3, 3]
    expect_equal(sqrt(vl), r1$coefficients$std.error)
    expect_equal(unname(rl$coefficients[3]), r1$coefficients$estimate)

    expect_equal(r1$coefficients$eff.obs, 954)
    expect_equal(954, sum(abs(r1$d$X)<=18))
    r1o <- capture.output(print(r1, digits=4))
    r2o <- capture.output(print(r2, digits=4))
    expect_equal(r1$coefficients$estimate, -1.1982581306)
    expect_equal(c(r1$coefficients$std.error, r2$coefficients$std.error),
                 c(0.6608206598, 0.6955768728))
    m1 <- paste0(" Sharp RD parameter   -1.198     ",
                 "0.6608        4.796     (-7.081, 4.684)")
    expect_equal(r1o[8], m1)
    expect_equal(r1o[10], "Onesided CIs:  (-Inf, 4.684), (-7.081, Inf)")
    expect_equal(r2o[10], "Onesided CIs:  (-Inf, 4.742), (-7.138, Inf)")
    expect_equal(r1o[16], "24 observations with missing values dropped")

    d <- NPRData(headst[!is.na(headst$mortHS), c("mortHS", "povrate")],
                 cutoff=0, "SRD")
    d <- PrelimVar(d, se.initial="Silverman")
    es <- function(kern, se.method) {
        NPRHonest(d, M=0.0076085544, kern=kern, sclass="H", se.method=se.method,
                  J=3, alpha=0.05, opt.criterion="MSE")
    }
    ff <- function(h, kern, se.method) {
        NPRHonest(d, M=0.0076085544, kern=kern, sclass="H", se.method=se.method,
                  J=3, alpha=0.05, h=h)
    }

    ## In this case the objective is not unimodal, but the modification still
    ## finds the minimum
    r <- es("uniform", "supplied.var")
    r1 <- es("uniform", "nn")
    ## Old algorithm
    ## expect_equal(unname(r$lower), -2.716800146662)
    ## expect_equal(unname(r1$estimate-r1$hl), -2.7294671445)
    ## Using SilvermanNN
    ## expect_equal(unname(r$lower), -2.6908821020)
    ## expect_equal(unname(r$h), 17.5696563721)
    ## expect_equal(unname(r1$estimate-r1$hl), -2.7106778586)
    expect_equal(r$coefficients$conf.low.onesided, -2.60568291)
    expect_equal(r$coefficients$bandwidth, 17.46128082)
    expect_equal(r1$coefficients$conf.low, -2.701699411)

    ## Just a sanity check
    expect_equal(r$coefficients$maximum.bias,
                 ff(r$coefficients$bandwidth,
                    "uniform", "supplied.var")$coefficients$maximum.bias)

    r <- es("triangular", "nn")
    expect_lt(abs(r$coefficients$bandwidth- 22.21108064), 5e-7)
    expect_lt(unname(r$coefficients$conf.high- 0.04129612), 1e-7)
    ## End replication

    ## Replicate 1511.06028v2, except we not longer allow se.initial=demeaned
    r <- RDHonest(voteshare ~ margin, data=lee08, kern="optimal",
                  M=0.2, opt.criterion="MSE", se.method="supplied.var")
    ## expect_equal(unname(r$lower), 2.2838100315)
    expect_equal(unname(r$coefficients$conf.low.onesided), 2.983141711)
    r2 <- RDHonest(voteshare ~ margin, data=lee08, kern="optimal", M=0.2,
                   opt.criterion="MSE", se.method="supplied.var",
                   alpha=r$coefficients$p.value)
    expect_equal(r2$coefficients$conf.low, 0)
    expect_equal(r$coefficients[1:4], r2$coefficients[1:4])

    r <- RDHonest(voteshare ~ margin, data=lee08, kern="optimal", M=0.04,
                  opt.criterion="OCI", se.method="supplied.var", beta=0.8)
    expect_equal(r$coefficients$conf.low.onesided, 2.761343298)
    r <- RDHonest(voteshare ~ margin, data=lee08, kern="optimal", M=0.04,
                  opt.criterion="FLCI", se.method="supplied.var")
    ro <- capture.output(print(r, digits=8))
    expect_equal(ro[8], paste0(" Sharp RD parameter 5.8864375  1.4472117",
                               "   0.80247208 (2.6646087, 9.1082664)"))
    ## Try nn
    r1 <- RDHonest(voteshare ~ margin, data=lee08, kern="optimal", M=0.04,
                   opt.criterion="FLCI", se.method="nn")
    expect_equal(as.numeric(r1$coefficients[2:6]),
                 c(5.886437523, 1.196647571, 0.8024720847, 3.099785428,
                   8.673089618))
    ## We no longer allow different bw below and above cutoff
    ## r <- RDHonest(voteshare ~ margin, data=lee08, kern="triangular",
    ##               M=0.04, opt.criterion="MSE", se.method="supplied.var",
    ##               se.initial="demeaned", bw.equal=FALSE, sclass="T")
    ## expect_equal(r$estimate, 5.9545948946)
    ## r1 <- capture.output(print(r, digits=6))
    ## expect_equal(r1[15],
    ##              "Bandwidth above cutoff: 10.6765")
    ## End replication

    ## Vignette replication Version: 0.1.2
    r <- RDHonest(voteshare ~ margin, data = lee08, kern = "uniform",
                  M = 0.1, h = 10, sclass = "T")
    expect_equal(r$coefficients$conf.low.onesided, 0.3162933187)
    r <- RDHonest(voteshare ~ margin, data = lee08, kern = "uniform",
                  M = 0.1, h = 10, sclass = "H")
    expect_equal(r$coefficients$conf.low.onesided, 2.3747626508)

    ## Again careful, not unimodal
    ## r1 <- RDOptBW(voteshare ~ margin, data = lee08, kern = "uniform",
    ##              M = 0.1, opt.criterion = "MSE", sclass = "T",
    ##              se.initial="SilvermanNN")
    ## r2 <- RDHonest(voteshare ~ margin, data = lee08, kern = "uniform",
    ##               M = 0.1, opt.criterion = "MSE", sclass = "H",
    ##               se.initial="SilvermanNN")
    ## r3 <- RDOptBW(voteshare ~ margin, data = lee08, kern = "uniform",
    ##              M = 0.1, opt.criterion = "MSE", sclass = "T",
    ##              se.initial="Silverman")
    ## Old optimization
    ## expect_equal(r1$h, 4.9367547102)
    ## expect_equal(unname(r2$lower), 2.3916108168)
    ## expect_equal(r3$h, 5.0590753991)

    ## Decrease M, these results are not true minima...
    d <- PrelimVar(NPRData(lee08, cutoff=0, "SRD"), se.initial="Silverman")
    r1 <- NPRHonest(d, M=0.01, kern="uniform", opt.criterion="MSE", beta=0.8,
                    sclass="T")$coefficients
    r2 <- NPRHonest(d, M=0.01, kern="uniform", opt.criterion="MSE", beta=0.8,
                    sclass="H")$coefficients
    r3 <- OptBW(d, kern = "uniform", M = 0.1, opt.criterion = "MSE",
                sclass = "T")
    expect_equal(unname(r1$bandwidth), 12.85186708)
    expect_equal(unname(r2$conf.low.onesided), 6.056860266)
    expect_equal(unname(r3), 5.086645484)

    ## Missing values
    expect_error(RDHonest(mortHS ~ povrate, data=headst, kern="uniform", h=12,
                          na.action="na.fail"))
    r1 <- RDHonest(mortHS ~ povrate, data=headst, kern="uniform",
                   na.action="na.omit")
    r1 <- capture.output(print(r1, digits=6))
    expect_equal(r1[c(11, 12, 16)],
                 c("Bandwidth: 3.98048, Kernel: uniform",
                   "Number of effective observations:     239",
                   "24 observations with missing values dropped"))
    r2 <- RDHonest(mortHS ~ povrate, data=headst, kern="epanechnikov",
                   point.inference=TRUE, cutoff=4)
    expect_equal(capture.output(print(r2, digits=6))[c(11, 12, 16)],
                 c("Bandwidth: 9.42625, Kernel: epanechnikov",
                   "Number of effective observations: 373.642",
                   "24 observations with missing values dropped"))
})

test_that("BME CIs match paper", {
    r1 <- RDHonestBME(log(earnings)~yearat14, data=cghs, cutoff=1947, h=6,
                      order=1)$coefficients
    expect_equal(as.numeric(r1[c("estimate", "std.error", "conf.low",
                                 "conf.high", "eff.obs")]),
                 c(0.0212923111, 0.03272404599, -0.13218978736, 0.17499392431,
                   20883))

    r2 <- RDHonestBME(log(earnings)~yearat14, cghs,
                      regformula="y~I(x>=0)+x+I(x^2)+I(x^3)+I(x^4)",
                      cutoff=1947)$coefficients
    ## Slightly different numbers with bug fixed
    expect_equal(as.numeric(r2[c("estimate", "std.error", "conf.low",
                                 "conf.high", "eff.obs")]),
                 c(0.05481510635, 0.02975117527, -0.2473386327, 0.3548003576,
                   73954))
    r3 <- RDHonestBME(log(cghs$earnings)~yearat14, data=cghs, h=3,
                      cutoff=1947, order=0, regformula="y~I(x>=0)")
    r4 <- RDHonestBME(log(cghs$earnings)~yearat14, data=cghs, h=3, order=0,
                      cutoff=1947)
    expect_equal(r3$coefficients, r4$coefficients)
    r5 <- RDHonestBME(duration~age, data=rebp, subset = (period==1 & female==0),
                      order=1, h=Inf, cutoff=50)
    r6 <- RDHonestBME(duration~age, data=rebp, subset = (period==1 & female==0),
                      order=3, h=1, cutoff=50, alpha=0.1)
    r5 <- capture.output(print(r5, digits=5))
    est <- paste0(" Sharp RD parameter   14.798      ",
                  "2.234       24.315   (-31.823, 60.724)")
    expect_equal(r5[8], est)
    expect_equal(r5[10], "Onesided CIs:  (-Inf, 57.249), (-28.236, Inf)")
    expect_equal(as.numeric(r6$coefficients[c(2:6, 10)]),
                 c(12.206023620, 8.878539149, 17.030704793, -24.727726605,
                   46.856041887, 3030))
})

test_that("Optimizing bw", {
    xprobs <- c(rep(.5/5, 5), rep(.5/4, 4))
    xsupp <- sort(c(-(1:5)/5, (1:4)/4))
    set.seed(42)
    x <- sample(xsupp, 100, prob=xprobs, replace=TRUE)
    d <- data.frame(y=rnorm(100, sd=1), x=x)

    r <- RDHonest(y~x, data=d, cutoff=0, M=40, kern="uniform",
                  opt.criterion="FLCI")
    expect_equal(r$coefficients$bandwidth, 1/2)
    r <- RDHonest(y~x, data=d, cutoff=0, M=2, kern="uniform",
                  opt.criterion="FLCI")
    expect_equal(r$coefficients$bandwidth, 0.8)

    xprobs <- c(rep(.5/4, 4), rep(.5/4, 4))
    xsupp <- sort(c(-(1:4)/4, (1:4)/4))
    set.seed(42)
    x <- sample(xsupp, 100, prob=xprobs, replace=TRUE)
    d <- data.frame(y=rnorm(100, sd=1), x=x)
    r1 <- RDHonest(y~x, data=d, cutoff=0, M=40, kern="triangular",
                   opt.criterion="FLCI")
    r2 <- RDHonest(y~x, data=d, cutoff=0, M=40, kern="triangular", h=0.51)
    expect_equal(r1$coefficients$conf.low, r2$coefficients$conf.low)
})

test_that("Supplied variance", {
    r <- RDHonest(voteshare ~ margin, data=lee08, M=2*0.00002, h=3)
    r2 <- RDHonest(I(10*voteshare) ~ margin, data=lee08, M=2*0.00002, h=3,
                   sigmaY2=r$data$sigma2, se.method="supplied.var")
    expect_equal(r$coefficients$std.error, r2$coefficients$std.error)
    r <- RDHonest(voteshare ~ margin, data=lee08, M=2*0.002, h=4,
                  point.inference=TRUE, cutoff=2)
    r2 <- RDHonest(I(10*voteshare) ~ margin, data=lee08, M=2*0.002, h=4,
                   sigmaY2=r$data$sigma2, se.method="supplied.var",
                   point.inference=TRUE, cutoff=2)
    expect_equal(r$coefficients$std.error, r2$coefficients$std.error)

    r <- RDHonest(log(cn)~retired | elig_year, data=rcp[1:100, ], cutoff=0,
                  M=c(0.002, 0.005), T0=0, h=7, kern="uniform")
    ## Data not in order
    r2 <- RDHonest(r$data$Y[, 1]~r$data$Y[, 2] | r$data$X, data=rcp[1:100, ],
                   cutoff=0, M=c(0.002, 0.005), T0=0, h=7,
                   sigmaY2=r$data$sigma2, kern="uniform")
    expect_equal(r$coefficients$std.error, r2$coefficients$std.error)


})


test_that("Adaptation bounds", {
    ## Replicate some analysis form ecta paper
    r <- RDHonest(voteshare ~ margin, data=lee08, M=2*0.00002, h=2)
    b1 <- RDTEfficiencyBound(r, opt.criterion="OCI", beta=0.8)
    r <- RDHonest(voteshare ~ margin, data=lee08, M=2*0.005, h=3)
    b2 <- RDTEfficiencyBound(r, 2*0.005, opt.criterion="FLCI")
    expect_equal(unname(c(b1, b2)), c(0.9956901, 0.95870418))
})


## rf <- RDHonest(log(cn)~retired | elig_year, data=rcp, cutoff=0)
## rs <- RDHonest(log(cn)~elig_year, data=rcp, cutoff=0)
## ro <- RDHonest(log(cn)~elig_year, data=rcp, cutoff=0, kern="optimal")
## rp <- RDHonest(log(cn)~elig_year, data=rcp, cutoff=1, point.inference=TRUE)
