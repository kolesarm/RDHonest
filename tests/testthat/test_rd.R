context("Test RD")

test_that("Test class constructor sorting", {
    d0 <- FRDData(rcp[, c(6, 3, 2)], cutoff=3)
    d1 <- FRDData(rcp[sort(rcp$elig_year,
                           index.return=TRUE)$ix, c(6, 3, 2)], cutoff=3)
    expect_identical(d1, d0)

    d0 <- RDData(rcp[, c(6, 2)], cutoff=3)
    d1 <- RDData(rcp[sort(rcp$elig_year,
                           index.return=TRUE)$ix, c(6, 2)], cutoff=3)
    expect_identical(d1, d0)

    d0 <- LPPData(rcp[, c(6, 2)], point=-3)
    d1 <- LPPData(rcp[sort(rcp$elig_year,
                           index.return=TRUE)$ix, c(6, 2)], point=-3)
    expect_identical(d1, d0)

    ## Now test that I() works
    expect_equal(RDHonest(voteshare ~ margin, data=lee08,
                          M=0, h=2)$coefficients$estimate,
                 -RDHonest(voteshare ~ I(-margin), data=lee08,
                           M=0, h=2)$coefficients$estimate)
    expect_equal(FRDHonest(cn~retired | elig_year, data=rcp, cutoff=0,
                           M=c(1, 0.1), h=3)$coefficients$estimate,
                 FRDHonest(cn~retired | I(2*elig_year), data=rcp, cutoff=0,
                           M=c(1, 0.1), h=6)$coefficients$estimate)
})

test_that("IK bandwidth calculations", {
    ## Test IK bandwidth in Lee data, IK Table 1
    d <- RDData(lee08, cutoff=0)

    dig <- getOption("digits")
    options(digits=8)
    r1 <- capture.output(IKBW.fit(d, verbose=TRUE))
    expect_equal(r1[c(2, 4, 5, 6, 8)],
                 c(" h1:  14.44507 ", " f(0):  0.0089622411 ",
                   " sigma^2_{+}(0):  12.024424 ^2",
                   " sigma^2_{+}(0): 10.472064 ^2 ",
                   " h_{2, +}: 60.513312 h_{2, -}: 60.993358 "))
    options(digits=dig)

    expect_equal(IKBW.fit(d), 29.3872649956)
    expect_equal(IKBW.fit(d, kern=function(u) pmax(1-abs(u), 0)), 29.3872649956)

    r <- NPRreg.fit(d, h=IKBW.fit(d, kern="uniform"), kern="uniform")
    expect_equal(r$estimate, 8.0770003749)
    d <- NPRPrelimVar.fit(RDData(lee08, cutoff=0), se.initial="EHW")
    expect_equal(sqrt(mean(d$sigma2p)), 12.58183131)
    expect_equal(sqrt(mean(d$sigma2m)), 10.79067278)
})

test_that("Plots", {
    if (requireNamespace("ggplot2", quietly = TRUE)) {
        expect_silent(invisible(plot_RDscatter(
            earnings~yearat14, data=cghs,
                   cutoff=1947, avg=Inf, propdotsize=TRUE)))
        expect_silent(invisible(plot_RDscatter(
            voteshare~margin, data=lee08, subset=abs(lee08$margin)<=50,
            avg=50, propdotsize=FALSE,
            xlab="margin", ylab="effect")))
        }
})


test_that("Honest inference in Lee and LM data",  {

    ## Replicate 1606.01200v2, except we no longer provide SilvermanNN
    r1 <- RDHonest(mortHS ~ povrate60, data=headst,
                  kern="uniform", h=18, M=0.1, se.method="EHW")
    r2 <- RDHonest(mortHS ~ povrate60, data=headst,
                  kern="uniform", h=18, M=0.1, se.method="nn")
    r1o <- capture.output(print(r1, digits=4))
    r2o <- capture.output(print(r2, digits=4))
    expect_equal(r1$coefficients$estimate, -1.1982581306)
    expect_equal(c(r1$coefficients$std.error, r2$coefficients$std.error),
                 c(0.6608206598, 0.6955768728))
    m1 <- paste0(" Sharp RD Parameter   -1.198     ",
                 "0.6608        4.796     (-7.081, 4.684)")
    expect_equal(r1o[8], m1)
    expect_equal(r1o[10], "Onesided CIs:  (-Inf, 4.684), (-7.081, Inf)")
    expect_equal(r2o[10], "Onesided CIs:  (-Inf, 4.742), (-7.138, Inf)")

    d <- RDData(headst[!is.na(headst$mortHS), c("mortHS", "povrate60")],
                cutoff=0)
    es <- function(kern, se.method) {
        NPRHonest.fit(d, M=0.0076085544, kern=kern, sclass="H",
                      se.method=se.method, J=3, alpha=0.05, opt.criterion="MSE",
                      se.initial="Silverman")
    }
    ff <- function(h, kern, se.method) {
        NPRHonest.fit(d, M=0.0076085544, kern=kern, sclass="H",
                      se.method=se.method, J=3, alpha=0.05, h=h,
                      se.initial="Silverman")
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
                  M=0.2, opt.criterion="MSE", se.method="supplied.var",
                  se.initial="EHW")
    ## expect_equal(unname(r$lower), 2.2838100315)
    expect_equal(unname(r$coefficients$conf.low.onesided), 2.983141711)
    r <- RDHonest(voteshare ~ margin, data=lee08, kern="optimal",
             M=0.04, opt.criterion="OCI", se.method="supplied.var",
             se.initial="EHW", beta=0.8)
    expect_equal(r$coefficients$conf.low.onesided, 2.761343298)
    r <- RDHonest(voteshare ~ margin, data=lee08, kern="optimal", M=0.04,
                  opt.criterion="FLCI", se.method="supplied.var",
                  se.initial="EHW")
    ro <- capture.output(print(r, digits=8))
    expect_equal(ro[9], paste0(" Sharp RD parameter 5.8864375  1.4472117",
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
    r1 <- RDHonest(voteshare ~ margin, data = lee08, kern = "uniform",
                 M = 0.01, opt.criterion = "MSE", sclass = "T",
                 se.initial="Silverman")$coefficients
    r2 <- RDHonest(voteshare ~ margin, data = lee08, kern = "uniform",
                  M = 0.01, opt.criterion = "MSE", sclass = "H",
                  se.initial="Silverman")$coefficients
    r3 <- NPROptBW.fit(RDData(lee08, cutoff=0), kern = "uniform", M = 0.1,
                       opt.criterion = "MSE", sclass = "T",
                       se.initial="Silverman")
    expect_equal(unname(r1$bandwidth), 12.85186708)
    expect_equal(unname(r2$conf.low.onesided), 6.056860266)
    expect_equal(unname(r3$h), 5.086645484)
})

test_that("BME CIs match paper", {
    r1 <- RDHonestBME(log(earnings)~yearat14, data=cghs,
                     cutoff=1947, h=6, order=1)
    expect_equal(as.numeric(r1$coefficients[
                                   c("estimate", "std.error", "conf.low",
                                     "conf.high", "eff.obs")]),
                 c(0.0212923111, 0.03272404599, -0.13218978736, 0.17499392431,
                   20883))

    r2 <- RDHonestBME(log(earnings)~yearat14, cghs, cutoff=1947,
                      regformula="y~I(x>=0)+x+I(x^2)+I(x^3)+I(x^4)")
    ## expect_equal(r1$CI, c(-0.23749230603, 0.34429708773))
    ## Slightly different numbers with bug fixed
    expect_equal(as.numeric(r2$coefficients[
                                   c("estimate", "std.error", "conf.low",
                                     "conf.high", "eff.obs")]),
                 c(0.05481510635, 0.02975117527, -0.2473386327, 0.3548003576,
                   73954))
    r3 <- RDHonestBME(log(cghs$earnings)~yearat14, data=cghs, h=3,
                      cutoff=1947, order=0, regformula="y~I(x>=0)")
    r4 <- RDHonestBME(log(cghs$earnings)~yearat14, data=cghs, h=3, order=0,
                      cutoff=1947)
    expect_equal(r3$coefficients, r4$coefficients)
    r5 <- RDHonestBME(duration~age, data=rebp, subset=(period==1 & female==0),
                      order=1, h=Inf, cutoff=50)
    r6 <- RDHonestBME(duration~age, data=rebp, subset=(period==1 & female==0),
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

    r <- RDHonest(y~x, data=d, cutoff=0, M=60, kern="uniform",
                  opt.criterion="FLCI", order=2)
    expect_equal(r$coefficients$bandwidth, 3/4)

    r <- RDHonest(y~x, data=d, cutoff=0, M=0.4, kern="uniform",
                  opt.criterion="FLCI", order=2)
    expect_equal(r$coefficients$bandwidth, 1)

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

test_that("Adaptation bounds", {
    ## Replicate some analysis form ecta paper
    r <- RDHonest(voteshare ~ margin, data=lee08, M=2*0.00002, h=2)
    b1 <- RDTEfficiencyBound(r, opt.criterion="OCI", beta=0.8)
    r <- RDHonest(voteshare ~ margin, data=lee08, M=2*0.005, h=3)
    b2 <- RDTEfficiencyBound(r, 2*0.005, opt.criterion="FLCI")
    expect_equal(unname(c(b1, b2)), c(0.9956901, 0.95870418))
})
