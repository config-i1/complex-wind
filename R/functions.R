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

#### cARIMA via clm, d = 1 ####
# clm's default starting values can end in a poor local optimum (a nested bigger model with a
# lower likelihood than the smaller one). So fit from the default start and from each supplied
# vector of starting values (complex: intercept, AR, MA), and keep the fit with the highest logLik.
fitCARIMA <- function(y, p, q, starts = list(), distribution = "dcnorm", drift = TRUE) {
  data <- data.frame(y = y)
  # drift = FALSE: no intercept (clm warns that this is "naughty"; it is intended here)
  formula <- if (drift) y ~ 1 else y ~ -1
  fit <- function(...) tryCatch(suppressWarnings(clm(formula, data, orders = c(p, 1, q),
                                                     distribution = distribution, ...)),
                                error = function(e) NULL)
  fits <- c(list(fit()), lapply(starts, function(B) fit(B = B)))
  fits <- fits[!sapply(fits, is.null)]
  if (length(fits) == 0) stop("all clm fits failed")
  fits[[which.max(sapply(fits, function(f) as.numeric(logLik(f))))]]
}
# Starting values for cARIMA(p, 1, q) from a nested model's coefficients (intercept if any, AR, MA),
# padding the new AR / MA terms with zeros
padStart <- function(coefs, pFrom, qFrom, p, q, intercept = TRUE) {
  k <- as.integer(intercept)
  c(coefs[seq_len(k)], coefs[k + seq_len(pFrom)], rep(0 + 0i, p - pFrom),
    coefs[k + pFrom + seq_len(qFrom)], rep(0 + 0i, q - qFrom))
}
forecastCARIMA <- function(fit, h) {
  as.vector(predict(fit, newdata = matrix(NA, h, 1))$mean)
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
