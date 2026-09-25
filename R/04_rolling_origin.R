# Rolling-origin point forecasts, h = 1..36, orders fixed from R/03_select_orders.R.
# Parameters are re-estimated at every origin on a 4-week window.
# Output: results/forecasts/<model>.rds, a complex matrix origin x horizon with the forecasts of
# Z = v * exp(i * theta), and results/forecasts/actuals.rds (same layout).
# Usage: Rscript R/04_rolling_origin.R [model ...]   (no arguments = all models)
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

# Each model: function(s) returning h complex forecasts of Z from the window series s
carimaModel <- function(name, series) {
  function(s) {
    o <- orders[[name]]
    # Starting values: the default and the estimates from the selection window (with the shape)
    start <- c(o$coef, if (o$distribution == "dcgnorm") o$shape)
    fc <- forecastCARIMA(fitCARIMA(s[[series]], o$p, o$q, list(start), o$distribution), h)
    if (series == "logZ") exp(fc) else fc
  }
}
realModel <- function(pkg, series1, series2, combine) {
  fc <- if (pkg == "arima") forecastArima else forecastMsarima
  function(s) {
    combine(fc(s[[series1]], orders[[paste0(pkg, "_", series1)]], h),
            fc(s[[series2]], orders[[paste0(pkg, "_", series2)]], h))
  }
}
models <- list(
  persistence = function(s) rep(tail(s$z, 1), h),
  carima_logZ = carimaModel("carima_logZ", "logZ"),
  carima_logZ_cgn = carimaModel("carima_logZ_cgn", "logZ"),
  carima_z = carimaModel("carima_z", "z"),
  carima_z_cgn = carimaModel("carima_z_cgn", "z")
)
for (pkg in c("arima", "msarima")) {
  models[[paste0(pkg, "_speed_angle")]] <- realModel(pkg, "speed", "angle", speedAngleZ)
  models[[paste0(pkg, "_logSpeed_angle")]] <- realModel(pkg, "logSpeed", "angle",
                                                        function(v, a) speedAngleZ(exp(v), a))
  models[[paste0(pkg, "_re_im")]] <- realModel(pkg, "re", "im",
                                               function(r, i) complex(real = r, imaginary = i))
}

run <- commandArgs(trailingOnly = TRUE)
if (length(run) == 0) run <- names(models)
stopifnot(all(run %in% names(models)))

startTime <- Sys.time()
results <- foreach(origin = setup$origins) %dopar% {
  s <- windowSeries(wind, (origin - setup$windowLength + 1):origin)
  sapply(run, function(m) tryCatch(models[[m]](s), error = function(e) rep(NA_complex_, h)))
}
for (m in run) {
  forecasts <- t(sapply(results, function(x) x[, m]))
  rownames(forecasts) <- setup$origins
  saveRDS(forecasts, paste0("results/forecasts/", m, ".rds"))
}
actuals <- t(sapply(setup$origins, function(origin) wind$z[origin + 1:h]))
rownames(actuals) <- setup$origins
saveRDS(actuals, "results/forecasts/actuals.rds")

cat("models:", paste(run, collapse = ", "), "\norigins:", length(setup$origins),
    " failed forecasts:", sum(sapply(results, function(x) sum(is.na(x[1, ])))),
    " time:", format(Sys.time() - startTime), "\n")
