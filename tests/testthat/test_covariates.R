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

test_that("Test collinear covariates", {

    df <- headst[complete.cases(headst), ]
    covlist <- "urban+black"
    fh1 <- as.formula(paste0("mortHS ~ povrate|", covlist))
    expect_message(r0 <- RDHonest(fh1, data=df), "Using Armstrong")
    Z <- cbind(df$urban, df$black, df$urban)
    expect_message(expect_message(r1 <- RDHonest(df$mortHS~df$povrate|Z,
                                                 M=r0$coefficients$M)))
    expect_lt(max(abs(r1$coefficients[2:13]-r0$coefficients[2:13])), 1e-15)
    Z <- cbind(df$urban, df$urban, df$black, df$black, df$urban+df$black)
    suppressMessages(expect_message(r2 <- RDHonest(df$mortHS~df$povrate|Z)))
    expect_lt(max(abs(r2$coefficients[2:13]-r0$coefficients[2:13])), 1e-15)
    expect_message(expect_message(r <- RDHonest(df$mortHS~df$povrate|Z, h=0.05,
                                                kern="uniform")), "large")

    ## Fuzzy
    df <- rcp[1:1000, ]
    Z <- cbind(df$food, df$food+df$family_size, df$family_size)
    expect_message(r0 <- RDHonest(log(c)|retired ~ elig_year|Z[, 1:2], data=df,
                                  weights=survey_year), "Using Armstrong")
    suppressMessages(expect_message(r1 <- RDHonest(log(c)|retired ~ elig_year|Z,
                                                   data=df,
                                                   weights=survey_year)))
    expect_lt(max(abs(r1$coefficients[2:13]-r0$coefficients[2:13])), 1e-15)
})
