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
    dd <- c(r0$d$coefs_on_covariates, r0$coefficients$estimate)-
        c(m3$coefficients[-idx], m3$coefficients[3])
    expect_lt(max(abs(dd)), 1e-6)

    ## TODO: return lm object, so we can do normal inference.

    if(FALSE) {
    z <- lm.wfit(fh2, data=df, weights = pmax(1 - abs(povrate/h), 0))
    z$na.action <- attr(mf, "na.action")
    z$contrasts <- attr(x, "contrasts")
    z$xlevels <- .getXlevels(mt, mf)
    z$call <- cl
    z$terms <- mt
    z$model <- mf

        fit <- lm(cbind(y1, y2, y3) ~ x1 + x2 + x3, data = dat)


    }




    ## TODO: pass function as kern

    ## TODO: test sort w/ covariates

})
