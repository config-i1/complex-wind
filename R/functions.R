# Shared functions for the experiment scripts.
suppressMessages({
  library(complex)
  library(forecast)
  library(smooth)
})

# Wrap an angle to (-pi, pi]
wrapAngle <- function(x) atan2(sin(x), cos(x))

# Unwrap an angle series: cumulative sum of shortest-arc changes, starting at the first value.
# Differences of the result are the true angle changes (no +-2pi jumps).
unwrapAngle <- function(theta) cumsum(c(theta[1], wrapAngle(diff(theta))))

# Series used by the models, from a window of the wind data
windowSeries <- function(wind, idx) {
  theta <- unwrapAngle(wind$theta[idx])
  list(speed = wind$speed[idx],
       logSpeed = log(wind$speed[idx]),
       angle = theta,
       re = Re(wind$z[idx]),
       im = Im(wind$z[idx]),
       z = wind$z[idx],
       logZ = complex(real = log(wind$speed[idx]), imaginary = theta))
}

#### cARIMA via clm ####
# Type of the constant: "none", "real", "imaginary" or "both". Earlier results store the drift as
# TRUE/FALSE (TRUE = complex constant), also as strings when read from csv.
constantType <- function(drift) {
  if (is.logical(drift) || drift %in% c("TRUE", "FALSE")) {
    if (as.logical(drift)) "both" else "none"
  } else {
    match.arg(as.character(drift), c("none", "real", "imaginary", "both"))
  }
}
# clm's default starting values can end in a poor local optimum (a nested bigger model with a
# lower likelihood than the smaller one). So fit from the default start and from each supplied
# vector of starting values (complex: intercept, AR, MA), and keep the fit with the highest logLik.
# hr = TRUE also starts from the Hannan-Rissanen estimates (clm(fast = TRUE)); used in Experiment 2.
fitCARIMA <- function(y, p, q, starts = list(), distribution = "dcnorm", drift = TRUE, d = 1, hr = FALSE) {
  data <- data.frame(y = y)
  # drift: the type of the constant (see constantType()). "none" = no intercept (clm warns that
  # this is "naughty"; it is intended here); "real" / "imaginary" = that part of it only
  formula <- switch(constantType(drift),
                    none = y ~ -1, both = y ~ 1, real = y ~ 1 - Im(1), imaginary = y ~ 1 - Re(1))
  fit <- function(...) tryCatch(suppressWarnings(clm(formula, data, orders = c(p, d, q),
                                                     distribution = distribution, ...)),
                                error = function(e) NULL)
  fits <- c(list(fit()), if (hr) list(fit(fast = TRUE)), lapply(starts, function(B) fit(B = B)))
  fits <- fits[!sapply(fits, is.null)]
  if (length(fits) == 0) stop("all clm fits failed")
  fits[[which.max(sapply(fits, function(f) as.numeric(logLik(f))))]]
}
# Starting values for cARIMA(p, 1, q) from a nested model's coefficients (intercept if any, AR, MA),
# padding the new AR / MA terms with zeros
padStart <- function(coefs, pFrom, qFrom, p, q, intercept = TRUE) {
  k <- as.integer(constantType(intercept) != "none")
  c(coefs[seq_len(k)], coefs[k + seq_len(pFrom)], rep(0 + 0i, p - pFrom),
    coefs[k + pFrom + seq_len(qFrom)], rep(0 + 0i, q - qFrom))
}
forecastCARIMA <- function(fit, h) {
  as.vector(predict(fit, newdata = matrix(NA, h, 1))$mean)
}
# complex::auto.carima chooses d (0 or 1), the constant (none, real, imaginary or both) and
# p, q (0..3) with its stepwise search.
# The result has the layout of the other cARIMA entries of results/orders.rds, plus d and the ICs.
selectAutoCarima <- function(y, distribution) {
  fit <- auto.carima(y, orders = list(ar = 3, i = 1, ma = 3), constant = "auto",
                     distribution = distribution)
  o <- fit$other$orders   # c(p, d, q); orders() does not work on carima objects
  list(p = o[1], d = o[2], q = o[3], drift = fit$constant, constant = fit$constant,
       distribution = distribution, coef = coef(fit),
       shape = if (is.null(fit$other$shape)) 1 else fit$other$shape, ICs = fit$ICs)
}

#### forecast::auto.arima / Arima, d = 1 ####
selectArima <- function(y) {
  fit <- auto.arima(y, d = 1, seasonal = FALSE)
  list(p = fit$arma[1], q = fit$arma[2], drift = "drift" %in% names(coef(fit)))
}
forecastArima <- function(y, spec, h) {
  fit <- Arima(y, order = c(spec$p, 1, spec$q), include.drift = spec$drift)
  as.numeric(forecast::forecast(fit, h = h)$mean)
}

#### smooth::auto.msarima / msarima, d = 1 ####
# auto.msarima also selects the order of integration. To keep d = 1, select AR/MA orders
# (and the constant) on the differenced series, then refit with i = 1 on the levels.
selectMsarima <- function(y) {
  fit <- auto.msarima(diff(y), orders = list(ar = 3, i = 0, ma = 3), lags = 1)
  o <- orders(fit)
  list(p = o$ar, q = o$ma, drift = !is.null(fit$constant))
}
forecastMsarima <- function(y, spec, h) {
  fit <- msarima(y, orders = list(ar = spec$p, i = 1, ma = spec$q), lags = 1,
                 constant = spec$drift)
  as.numeric(forecast(fit, h = h)$mean)
}
