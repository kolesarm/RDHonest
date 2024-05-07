test_that("Test inputs", {
    expect_error(r2 <- RDHonest(log(cn)~retired|elig_year, data=rcp, M=1, T0=0))
    expect_error(r2 <- RDHonest(log(cn)~elig_year, data=rcp, M=c(1, 1)))
    expect_error(r2 <- RDHonest(log(cn)~elig_year, data=rcp, kern="Unif"))
    expect_error(RDHonest(log(cn)|retired~elig_year, data=rcp, M=1, T0=0,
                          kern="optimal"))
    expect_error(pp <- RDHonest(voteshare~margin, data=lee08,
                                M=2, h=5, subset=I(margin>0))$coefficients)
    expect_error(RDHonest(c~elig_year, data=rcp, clusterid=seq_along(rr$c),
                          point.inference=TRUE))
    ## Insufficient unique values of the running variable to compute rule of
    ## thumb for M.
    expect_error(r2 <- RDHonest(log(cn)~elig_year,
                                data=rcp[abs(rcp$elig_year)<=3, ]))
    expect_error(RDHonest(cn ~ retired | elig_year, data=rcp))
})
