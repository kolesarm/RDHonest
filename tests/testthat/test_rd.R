context("Test RD")

test_that("IK bandwidth calculations", {
    ## Test IK bandwidth in Lee data, IK Table 1
    d <- RDData(lee08, cutoff=0)
    expect_equal(IKBW.fit(d), 29.3872649956)

    r <- RDLPreg(d, hp=IKBW.fit(d, kern="uniform"), kern="uniform")
    expect_equal(r$estimate, 8.0770003749)
    d <- RDprelimVar(RDData(lee08, cutoff=0), se.initial="IK")
    expect_equal(sqrt(mean(d$sigma2p)), 14.5030086128)
    expect_equal(sqrt(mean(d$sigma2m)), 12.4627517531)
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
                     se.initial="SilvermanOld")
    }
    r <- es("uniform", "supplied.var")
    expect_equal(r$hm, 18)
    r <- es("uniform", "nn")
    expect_equal(unname(r$estimate-r$hl), -2.7294671445)
    r <- es("triangular", "nn")
    expect_equal(r$hm, 22.8088245304)
    expect_equal(unname(r$estimate+r$hl), 0.0547661086)
    ## End replication

    ## Replicate 1511.06028v2
    r <- RDHonest(voteshare ~ margin, data=lee08, kern="optimal",
                  M=0.2, opt.criterion="MSE", se.method="supplied.var",
                  se.initial="IK")
    expect_equal(unname(r$lower), 2.2838100315)
    r <- RDHonest(voteshare ~ margin, data=lee08, kern="optimal",
             M=0.04, opt.criterion="OCI", se.method="supplied.var",
             se.initial="IK", beta=0.8)
    expect_equal(r$lower, c("supplied.var"=2.2802786871))
    r <- RDHonest(voteshare ~ margin, data=lee08, kern="optimal",
                  M=0.04, opt.criterion="FLCI", se.method="supplied.var",
                  se.initial="IK")
    expect_equal(r$estimate-r$hl, c("supplied.var"=2.3786060461))

    r <- RDHonest(voteshare ~ margin, data=lee08, kern="triangular",
                  M=0.04, opt.criterion="MSE", se.method="supplied.var",
                  se.initial="IK", bw.equal=FALSE, sclass="T")
    expect_equal(r$estimate, 5.9545948946)
    ## End replication

    ## Vignette replication Version: 0.1.2
    r <- RDHonest(voteshare ~ margin, data = lee08, kern = "uniform",
                  M = 0.1, hp = 10, sclass = "T")
    expect_equal(unname(r$lower), 0.3162933187)
    r <- RDHonest(voteshare ~ margin, data = lee08, kern = "uniform",
                  M = 0.1, hp = 10, sclass = "H")
    expect_equal(unname(r$lower), 2.3747626508)

    r <- RDOptBW(voteshare ~ margin, data = lee08, kern = "uniform",
                 M = 0.1, opt.criterion = "MSE", sclass = "T",
                 se.initial="SilvermanOld")
    expect_equal(r$hp, 4.9367547102)
    r <- RDHonest(voteshare ~ margin, data = lee08, kern = "uniform",
                  M = 0.1, opt.criterion = "MSE", sclass = "H",
                  se.initial="SilvermanOld")
    expect_equal(unname(r$lower), 2.3916108168)
    r <- RDOptBW(voteshare ~ margin, data = lee08, kern = "uniform",
                 M = 0.1, opt.criterion = "MSE", sclass = "T",
                 se.initial="Silverman")
    expect_equal(r$hp, 5.0590753991)
})
