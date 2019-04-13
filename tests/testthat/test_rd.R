context("Test RD")

test_that("IK bandwidth calculations", {
    ## Test IK bandwidth in Lee data, IK Table 1
    d <- RDData(lee08, cutoff=0)
    expect_equal(IKBW.fit(d), 29.3872649956)

    r <- RDLPreg(d, hp=IKBW.fit(d, kern="uniform"), kern="uniform")
    expect_equal(r$estimate, 8.0770003749)
    d <- RDPrelimVar(RDData(lee08, cutoff=0), se.initial="IKdemeaned")
    expect_equal(sqrt(mean(d$sigma2p)), 14.5030086128)
    expect_equal(sqrt(mean(d$sigma2m)), 12.4627517531)
})

test_that("Plots", {
    if (requireNamespace("ggplot2", quietly = TRUE)) {
        expect_silent(invisible(plot_RDscatter(
            RDData(data.frame(y=log(cghs$earnings),
                              x=cghs$yearat14),
                   cutoff=1947), avg=Inf, propdotsize=TRUE)))
        }
})


test_that("Honest inference in Lee and LM data",  {

    ## Replicate 1606.01200v2
    r <- RDHonest(mortHS ~ povrate60, data=headst,
                  kern="uniform", hp=18, M=0.1, se.method=c("EHW", "nn"))
    expect_equal(r$estimate, -1.1982581306)
    expect_equal(unname(r$sd), c(0.6608206598, 0.6955768728))

    d <- RDData(headst[!is.na(headst$mortHS), c("mortHS", "povrate60")],
                cutoff=0)
    es <- function(kern, se.method) {
        RDHonest.fit(d, M=0.0076085544, kern=kern, sclass="H",
                     se.method=se.method,
                     J=3, alpha=0.05, opt.criterion="MSE", bw.equal=TRUE,
                     se.initial="SilvermanNN")
    }
    ff <- function(h, kern, se.method) {
        RDHonest.fit(d, M=0.0076085544, kern=kern, sclass="H",
                     se.method=se.method,
                     J=3, alpha=0.05, hp=h, bw.equal=TRUE,
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
    expect_equal(unname(r$hp), 17.5696563721)
    expect_equal(unname(r1$estimate-r1$hl), -2.7106778586)

    ## Just a sanity check
    expect_equal(r$maxbias, ff(r$hp, "uniform", "supplied.var")$maxbias)

    r <- es("triangular", "nn")
    expect_equal(r$hm, 22.8088245304)
    expect_equal(unname(r$estimate+r$hl), 0.0547661086)
    ## End replication

    ## Replicate 1511.06028v2
    r <- RDHonest(voteshare ~ margin, data=lee08, kern="optimal",
                  M=0.2, opt.criterion="MSE", se.method="supplied.var",
                  se.initial="IKdemeaned")
    expect_equal(unname(r$lower), 2.2838100315)
    r <- RDHonest(voteshare ~ margin, data=lee08, kern="optimal",
             M=0.04, opt.criterion="OCI", se.method="supplied.var",
             se.initial="IKdemeaned", beta=0.8)
    expect_equal(r$lower, c("supplied.var"=2.2802786871))
    r <- RDHonest(voteshare ~ margin, data=lee08, kern="optimal",
                  M=0.04, opt.criterion="FLCI", se.method="supplied.var",
                  se.initial="IKdemeaned")
    expect_equal(r$estimate-r$hl, c("supplied.var"=2.3786060461))

    r <- RDHonest(voteshare ~ margin, data=lee08, kern="triangular",
                  M=0.04, opt.criterion="MSE", se.method="supplied.var",
                  se.initial="IKdemeaned", bw.equal=FALSE, sclass="T")
    expect_equal(r$estimate, 5.9545948946)
    ## End replication

    ## Vignette replication Version: 0.1.2
    r <- RDHonest(voteshare ~ margin, data = lee08, kern = "uniform",
                  M = 0.1, hp = 10, sclass = "T")
    expect_equal(unname(r$lower), 0.3162933187)
    r <- RDHonest(voteshare ~ margin, data = lee08, kern = "uniform",
                  M = 0.1, hp = 10, sclass = "H")
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
    ## expect_equal(r1$hp, 4.9367547102)
    ## expect_equal(unname(r2$lower), 2.3916108168)
    ## expect_equal(r3$hp, 5.0590753991)

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

    expect_equal(r1$hp, 12.6576125622)
    expect_equal(unname(r2$lower), 6.0484981004)
    expect_equal(r3$hp, 5.0866454840)

})

test_that("BME CIs match paper", {
    ## Test IK bandwidth in Lee data, IK Table 1
    r1 <- RDHonestBME(log(earnings)~yearat14, data=cghs,
                     cutoff=1947, hp=6, order=1)
    expect_equal(r1$CI, c(-0.13218978736, 0.17499392431))
    r1 <- RDHonestBME(log(earnings)~yearat14, cghs, cutoff=1947,
                      regformula="y~I(x>=0)+x+I(x^2)+I(x^3)+I(x^4)")

    expect_equal(r1$CI, c(-0.23749230603, 0.34429708773))
})

test_that("Optmizing bw", {
    xprobs <- c(rep(.5/5, 5), rep(.5/4, 4))
    xsupp <- sort(c(-(1:5)/5, (1:4)/4))
    set.seed(42)
    x <- sample(xsupp, 100, prob=xprobs, replace=TRUE)
    d <- data.frame(y=rnorm(100, sd=1), x=x)

    r <- RDHonest(y~x, data=d, cutoff=0, M=40, kern="uniform",
                  opt.criterion="FLCI")
    expect_equal(r$hp, 1/2)

    r <- RDHonest(y~x, data=d, cutoff=0, M=60, kern="uniform",
                  opt.criterion="FLCI", order=2)
    expect_equal(r$hp, 3/4)

    r <- RDHonest(y~x, data=d, cutoff=0, M=0.4, kern="uniform",
                  opt.criterion="FLCI", order=2)
    expect_equal(r$hp, 1)

    xprobs <- c(rep(.5/4, 4), rep(.5/4, 4))
    xsupp <- sort(c(-(1:4)/4, (1:4)/4))
    set.seed(42)
    x <- sample(xsupp, 100, prob=xprobs, replace=TRUE)
    d <- data.frame(y=rnorm(100, sd=1), x=x)
    r1 <- RDHonest(y~x, data=d, cutoff=0, M=40, kern="triangular",
                  opt.criterion="FLCI")
    r2 <- RDHonest(y~x, data=d, cutoff=0, M=40, kern="triangular", hp=0.51)
    expect_equal(r1$lower, r2$lower)


})
