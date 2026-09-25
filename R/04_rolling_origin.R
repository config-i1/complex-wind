# Rolling-origin point forecasts, h = 1..36, orders fixed from R/03_select_orders.R.
# Parameters are re-estimated at every origin on a 4-week window.
# Output: results/forecasts/forecasts.rds, a list (one element per origin) of complex matrices
# h x model holding the forecast of Z = v * exp(i * theta), plus the actuals.
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
source("R/functions.R")
library(doMC)
nWorkers <- 30
registerDoMC(nWorkers)

wind <- readRDS("data/processed/wind.rds")
setup <- readRDS("data/processed/origins.rds")
orders <- readRDS("results/orders.rds")
h <- setup$horizon

# Build Z from separate forecasts of speed and angle. A negative speed forecast is set to zero.
speedAngleZ <- function(speed, angle) complex(modulus = pmax(speed, 0), argument = angle)

forecastOrigin <- function(origin) {
  idx <- (origin - setup$windowLength + 1):origin
  s <- windowSeries(wind, idx)
  out <- list()

  # Persistence
  out$persistence <- rep(tail(s$z, 1), h)

  # cARIMA on log(z) = log(speed) + i angle; exp() returns to the original units
  # Starting values: the default and the estimates from the selection window
  o <- orders$carima_logZ
  out$carima_logZ <- tryCatch(exp(forecastCARIMA(fitCARIMA(s$logZ, o$p, o$q, list(o$coef)), h)),
                              error = function(e) rep(NA_complex_, h))

  # cARIMA on z
  o <- orders$carima_z
  out$carima_z <- tryCatch(forecastCARIMA(fitCARIMA(s$z, o$p, o$q, list(o$coef)), h),
                           error = function(e) rep(NA_complex_, h))

  # Independent real-valued models
  for (pkg in c("arima", "msarima")) {
    fc <- if (pkg == "arima") forecastArima else forecastMsarima
    f <- function(series) tryCatch(fc(s[[series]], orders[[paste0(pkg, "_", series)]], h),
                                   error = function(e) rep(NA_real_, h))
    angle <- f("angle")
    out[[paste0(pkg, "_speed_angle")]] <- speedAngleZ(f("speed"), angle)
    out[[paste0(pkg, "_logSpeed_angle")]] <- speedAngleZ(exp(f("logSpeed")), angle)
    out[[paste0(pkg, "_re_im")]] <- complex(real = f("re"), imaginary = f("im"))
  }

  list(origin = origin, forecast = do.call(cbind, out), actual = wind$z[origin + 1:h])
}

startTime <- Sys.time()
forecasts <- foreach(origin = setup$origins, .errorhandling = "pass") %dopar% forecastOrigin(origin)
saveRDS(forecasts, "results/forecasts/forecasts.rds")
cat("origins:", length(forecasts), " failed:", sum(sapply(forecasts, inherits, "error")),
    " time:", format(Sys.time() - startTime), "\n")
