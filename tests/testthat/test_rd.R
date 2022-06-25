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
    expect_equal(RDHonest(voteshare ~ margin, data=lee08, M=0, h=2)$estimate,
                 -RDHonest(voteshare ~ I(-margin),
                           data=lee08, M=0, h=2)$estimate)
    expect_equal(FRDHonest(cn~retired | elig_year, data=rcp, cutoff=0,
                           M=c(1, 0.1), h=3)$estimate,
                 FRDHonest(cn~retired | I(2*elig_year), data=rcp, cutoff=0,
                           M=c(1, 0.1), h=6)$estimate)
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
    d <- NPRPrelimVar.fit(RDData(lee08, cutoff=0), se.initial="demeaned")
    expect_equal(sqrt(mean(d$sigma2p)), 14.5030086128)
    expect_equal(sqrt(mean(d$sigma2m)), 12.4627517531)
})

test_that("Plots", {
    if (requireNamespace("ggplot2", quietly = TRUE)) {
        expect_silent(invisible(plot_RDscatter(
            RDData(data.frame(y=log(cghs$earnings),
                              x=cghs$yearat14),
                   cutoff=1947), avg=Inf, propdotsize=TRUE)))
        expect_silent(invisible(plot_RDscatter(
            RDData(lee08, cutoff=0), avg=50, propdotsize=FALSE,
            window=50, xlab="margin", ylab="effect")))
        }
})


test_that("Honest inference in Lee and LM data",  {

    ## Replicate 1606.01200v2
    r <- RDHonest(mortHS ~ povrate60, data=headst,
                  kern="uniform", h=18, M=0.1, se.method=c("EHW", "nn"))
    r0 <- capture.output(print(r, digits=4))
    expect_equal(r$estimate, -1.1982581306)
    expect_equal(unname(r$sd), c(0.6608206598, 0.6955768728))
    expect_equal(r0[13], "nn    (-7.138, 4.742), (-7.138, Inf), (-Inf, 4.742)")

    d <- RDData(headst[!is.na(headst$mortHS), c("mortHS", "povrate60")],
                cutoff=0)
    es <- function(kern, se.method) {
        NPRHonest.fit(d, M=0.0076085544, kern=kern, sclass="H",
                      se.method=se.method, J=3, alpha=0.05, opt.criterion="MSE",
                      se.initial="SilvermanNN")
    }
    ff <- function(h, kern, se.method) {
        NPRHonest.fit(d, M=0.0076085544, kern=kern, sclass="H",
                      se.method=se.method, J=3, alpha=0.05, h=h,
                      se.initial="SilvermanNN")
    }

    ## In this case the objective is not unimodal, but the modification still
    ## finds the minimum
    r <- es("uniform", "supplied.var")
    r1 <- es("uniform", "nn")
    ## Old algorithm
    ## expect_equal(unname(r$lower), -2.716800146662)
    ## expect_equal(unname(r1$estimate-r1$hl), -2.7294671445)
    ## New algorithm
    expect_equal(unname(r$lower), -2.6908821020)
    expect_equal(unname(r$h), 17.5696563721)
    expect_equal(unname(r1$estimate-r1$hl), -2.7106778586)

    ## Just a sanity check
    expect_equal(r$maxbias, ff(r$h, "uniform", "supplied.var")$maxbias)

    r <- es("triangular", "nn")
    expect_lt(abs(r$h- 22.80882408), 5e-7)
    expect_lt(unname(r$estimate+r$hl- 0.05476609), 1e-7)
    ## End replication

    ## Replicate 1511.06028v2
    r <- RDHonest(voteshare ~ margin, data=lee08, kern="optimal",
                  M=0.2, opt.criterion="MSE", se.method="supplied.var",
                  se.initial="demeaned")
    expect_equal(unname(r$lower), 2.2838100315)
    r <- RDHonest(voteshare ~ margin, data=lee08, kern="optimal",
             M=0.04, opt.criterion="OCI", se.method="supplied.var",
             se.initial="demeaned", beta=0.8)
    expect_equal(r$lower, c("supplied.var"=2.2802786871))
    r <- RDHonest(voteshare ~ margin, data=lee08, kern="optimal",
                  M=0.04, opt.criterion="FLCI", se.method="supplied.var",
                  se.initial="demeaned")
    expect_equal(r$estimate-r$hl, c("supplied.var"=2.3786060461))

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
    expect_equal(unname(r$lower), 0.3162933187)
    r <- RDHonest(voteshare ~ margin, data = lee08, kern = "uniform",
                  M = 0.1, h = 10, sclass = "H")
    expect_equal(unname(r$lower), 2.3747626508)

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
    r1 <- RDOptBW(voteshare ~ margin, data = lee08, kern = "uniform",
                 M = 0.01, opt.criterion = "MSE", sclass = "T",
                 se.initial="SilvermanNN")
    r2 <- RDHonest(voteshare ~ margin, data = lee08, kern = "uniform",
                  M = 0.01, opt.criterion = "MSE", sclass = "H",
                  se.initial="SilvermanNN")
    r3 <- RDOptBW(voteshare ~ margin, data = lee08, kern = "uniform",
                 M = 0.1, opt.criterion = "MSE", sclass = "T",
                 se.initial="Silverman")
    expect_equal(unname(r1$h), 12.6576125622)
    expect_equal(unname(r2$lower), 6.0484981004)
    expect_equal(unname(r3$h), 5.0866454840)
})

test_that("BME CIs match paper", {
    r1 <- RDHonestBME(log(earnings)~yearat14, data=cghs,
                     cutoff=1947, h=6, order=1)
    expect_equal(r1$CI, c(-0.13218978736, 0.17499392431))
    r2 <- capture.output(print(r1, digits=5))
    expect_equal(r2[7], "(-0.13219, 0.17499)")

    r1 <- RDHonestBME(log(earnings)~yearat14, cghs, cutoff=1947,
                      regformula="y~I(x>=0)+x+I(x^2)+I(x^3)+I(x^4)")
    expect_equal(r1$CI, c(-0.23749230603, 0.34429708773))
    r3 <- RDHonestBME(log(cghs$earnings)~yearat14, data=cghs, h=3,
                      cutoff=1947, order=0, regformula="y~I(x>=0)")
    r4 <- RDHonestBME(log(cghs$earnings)~yearat14, data=cghs, h=3, order=0,
                      cutoff=1947)
    expect_equal(r3$CI, r4$CI)
})

test_that("Optimizing bw", {
    xprobs <- c(rep(.5/5, 5), rep(.5/4, 4))
    xsupp <- sort(c(-(1:5)/5, (1:4)/4))
    set.seed(42)
    x <- sample(xsupp, 100, prob=xprobs, replace=TRUE)
    d <- data.frame(y=rnorm(100, sd=1), x=x)

    r <- RDHonest(y~x, data=d, cutoff=0, M=40, kern="uniform",
                  opt.criterion="FLCI")
    expect_equal(r$h, 1/2)

    r <- RDHonest(y~x, data=d, cutoff=0, M=60, kern="uniform",
                  opt.criterion="FLCI", order=2)
    expect_equal(r$h, 3/4)

    r <- RDHonest(y~x, data=d, cutoff=0, M=0.4, kern="uniform",
                  opt.criterion="FLCI", order=2)
    expect_equal(r$h, 1)

    xprobs <- c(rep(.5/4, 4), rep(.5/4, 4))
    xsupp <- sort(c(-(1:4)/4, (1:4)/4))
    set.seed(42)
    x <- sample(xsupp, 100, prob=xprobs, replace=TRUE)
    d <- data.frame(y=rnorm(100, sd=1), x=x)
    r1 <- RDHonest(y~x, data=d, cutoff=0, M=40, kern="triangular",
                  opt.criterion="FLCI")
    r2 <- RDHonest(y~x, data=d, cutoff=0, M=40, kern="triangular", h=0.51)
    expect_equal(r1$lower, r2$lower)
})

test_that("Adaptation bounds", {
    ## Replicate some analysis form ecta paper
    d <- RDData(lee08, cutoff=0)
    b1 <- RDTEfficiencyBound(d, 2*0.00002, opt.criterion="OCI", alpha=0.05,
                             beta=0.8)
    b2 <- RDTEfficiencyBound(d, 2*0.005, opt.criterion="FLCI", alpha=0.05)
    expect_equal(unname(c(b1, b2)), c(0.9956901, 0.95870418))
})
