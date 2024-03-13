test_that("Test covariates", {
    ## Sharp design, no M or h
    df <- headst[complete.cases(headst), ]
    covlist <- "urban*black+sch1417"
    fh1 <- as.formula(paste0("mortHS ~ povrate|", covlist))
    fh2 <- as.formula(paste0("mortHS~povrate*I(povrate>=0)+", covlist))
    expect_message(r0 <- RDHonest(fh1, data=df))
    ## Doing it manually
    expect_message(r1 <- RDHonest(mortHS ~ povrate, data=df))
    h <- r1$coefficients$bandwidth
    m2 <- lm(fh2, data=df, weights = pmax(1 - abs(povrate/h), 0))
    idx <- c(1:3, 7)
    df$Y <- df$mortHS -
        model.matrix(m2$model, data=df)[, -idx] %*% m2$coefficients[-idx]
    expect_message(r2 <- RDHonest(Y~povrate, data=df))
    h2 <- r2$coefficients$bandwidth
    m3 <- lm(fh2, data=df, weights = pmax(1 - abs(povrate/h2), 0))
    expect_lt(abs(r0$coefficients$bandwidth-h2), 1e-5)
    dd <- sort(r0$lm$coefficients)-sort(m3$coefficients)
    expect_lt(max(abs(dd)), 1e-9)
    ## with h /M
    r0 <- RDHonest(fh1, data=df, M=2)
    r2 <- RDHonest(Y~povrate, data=df, M=2)
    h2 <- r2$coefficients$bandwidth
    m3 <- lm(fh2, data=df, weights = pmax(1 - abs(povrate/h2), 0))
    expect_lt(abs(r0$coefficients$bandwidth-h2), 1e-5)
    dd <- sort(r0$lm$coefficients)-sort(m3$coefficients)
    expect_lt(max(abs(dd)), 1e-9)

    ## fuzzy
    df <- rcp[1:1000, ]
    expect_message(r0 <- RDHonest(log(c)|retired ~ elig_year|food, data=df,
                                  weights=survey_year))
    ## Doing it manually
    expect_message(r1 <- RDHonest(log(c)|retired ~ elig_year, data=df,
                                  weights=survey_year))
    h <- r1$coefficients$bandwidth
    m2 <- lm(cbind(log(c), retired) ~ elig_year*I(elig_year>0)+food,
             data=df, weights = pmax(1 - abs(elig_year/h), 0)*survey_year)
    Y <- cbind(log(df$c), df$retired) - df$food %o% m2$coefficients[4, ]
    expect_message(r2 <- RDHonest(Y[, 1]|Y[, 2] ~ elig_year, data=df,
                                  weights=survey_year))
    h2 <- r2$coefficients$bandwidth
    expect_lt(abs(r0$coefficients$bandwidth-h2), 1e-5)
    m3 <- lm(cbind(log(c), retired) ~ elig_year*I(elig_year>0)+food,
             data=df, weights = pmax(1 - abs(elig_year/h2), 0)*survey_year)
    dd <- sort(r0$lm$coefficients)-sort(m3$coefficients)
    expect_lt(max(abs(dd)), 1e-8)

    ## pass function as kern
    expect_message(r00 <- RDHonest(log(c)|retired ~ elig_year|food, data=df,
                                   weights=survey_year,
                                   h=r0$coefficients$bandwidth,
                                   kern=function(u) pmax(1-abs(u), 0)))
    expect_equal(r00$coefficients, r00$coefficients)
})
