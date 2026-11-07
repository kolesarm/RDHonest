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
    expect_lt(abs(r0$coefficients$bandwidth-h2), 1e-5)
    ## Compare the regression at r0's own bandwidth: the two bandwidths agree
    ## only to optimizer precision, which the check above covers
    h0 <- r0$coefficients$bandwidth
    m3 <- lm(fh2, data=df, weights = pmax(1 - abs(povrate/h0), 0))
    dd <- sort(r0$lm$coefficients)-sort(m3$coefficients)
    expect_lt(max(abs(dd)), 1e-9)
    ## with h /M
    r0 <- RDHonest(fh1, data=df, M=2)
    r2 <- RDHonest(Y~povrate, data=df, M=2)
    h2 <- r2$coefficients$bandwidth
    expect_lt(abs(r0$coefficients$bandwidth-h2), 1e-5)
    ## Compare the regression at r0's own bandwidth: the two bandwidths agree
    ## only to optimizer precision, which the check above covers
    h0 <- r0$coefficients$bandwidth
    m3 <- lm(fh2, data=df, weights = pmax(1 - abs(povrate/h0), 0))
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
    ## As above, compare the regression at r0's own bandwidth
    h0 <- r0$coefficients$bandwidth
    m3 <- lm(cbind(log(c), retired) ~ elig_year*I(elig_year>0)+food,
             data=df, weights = pmax(1 - abs(elig_year/h0), 0)*survey_year)
    dd <- (sort(r0$lm$coefficients)-sort(m3$coefficients)) /
        sort(r0$lm$coefficients)
    expect_lt(max(abs(dd)), 1e-5)

    ## pass function as kern
    expect_message(r00 <- RDHonest(log(c)|retired ~ elig_year|food, data=df,
                                   weights=survey_year,
                                   h=r0$coefficients$bandwidth,
                                   kern=function(u) pmax(1-abs(u), 0)))
    expect_message(r0h <- RDHonest(log(c)|retired ~ elig_year|food, data=df,
                                   weights=survey_year,
                                   h=r0$coefficients$bandwidth))
    ## Kernel constants for a user-supplied kernel may be computed numerically,
    ## so allow more than rounding error
    expect_equal(r00$coefficients[2:13], r0h$coefficients[2:13],
                 tolerance = tol_opt)
})

test_that("Collinear covariates are dropped correctly", {
    ## Sharp
    df <- headst[complete.cases(headst), ]
    expect_message(r0 <- RDHonest(mortHS ~ povrate | urban + black, data = df),
                   "Using Armstrong")
    M0 <- r0$coefficients$M
    h0 <- r0$coefficients$bandwidth

    ## Fix M and h at the reference values so that only the handling of
    ## collinear columns differs; results then agree up to linear-algebra
    ## rounding, and the default tolerance applies.
    Z1 <- cbind(df$urban, df$black, df$urban)
    expect_message(r1 <- RDHonest(df$mortHS ~ df$povrate | Z1, M = M0, h = h0))
    expect_equal(as.vector(r1$coefficients[2:13]),
                 as.vector(r0$coefficients[2:13]))

    Z2 <- cbind(df$urban, df$urban, df$black, df$black, df$urban + df$black)
    expect_message(r2 <- RDHonest(df$mortHS ~ df$povrate | Z2, M = M0, h = h0))
    expect_equal(as.vector(r2$coefficients[2:13]),
                 as.vector(r0$coefficients[2:13]))

    ## End to end: M and h re-selected. Selected values are determined only to
    ## about sqrt(.Machine$double.eps) relative, so allow for that.
    suppressMessages(r2s <- RDHonest(df$mortHS ~ df$povrate | Z2))
    expect_equal(as.vector(r2s$coefficients[2:13]),
                 as.vector(r0$coefficients[2:13]),
                 tolerance = tol_opt)

    expect_message(expect_message(RDHonest(df$mortHS ~ df$povrate | Z2,
                                           h = 0.05, kern = "uniform")),
                   "large")

    ## Fuzzy
    df <- rcp[1:1000, ]
    Z <- cbind(df$food, df$food + df$family_size, df$family_size)
    expect_message(r0 <- RDHonest(log(c) | retired ~ elig_year | Z[, 1:2],
                                  data = df, weights = survey_year),
                   "Using Armstrong")
    ## Reduced-form and first-stage smoothness constants
    M0 <- c(r0$coefficients$M.rf, r0$coefficients$M.fs)
    h0 <- r0$coefficients$bandwidth
    expect_message(r1 <- RDHonest(log(c) | retired ~ elig_year | Z, data = df,
                                  weights = survey_year, M = M0, h = h0))
    expect_equal(r1$coefficients[2:13], r0$coefficients[2:13])
})
