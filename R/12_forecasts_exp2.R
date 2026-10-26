# Experiment 2: point forecasts and predictive distributions (quantiles 0.01..0.99 of speed and
# direction) for h = 1..12 at the 300 test origins. Parameters re-estimated at every origin.
# Orders: grid cARIMA from cross-validation (results/exp2/orders_cv.rds); auto.carima models and
# real benchmarks as selected in Experiment 1 (results/tables/orders.csv).
# Usage: Rscript R/12_forecasts_exp2.R [model ...]   (no arguments = all models)
# Output: results/exp2/forecasts/<model>.rds: list(point = origin x h complex, speed and
# direction = origin x h x tau arrays), results/exp2/forecasts/actuals.rds
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
source("R/exp2_functions.R")
library(doMC)
registerDoMC(min(parallel::detectCores() - 2, 30))

wind <- readRDS("data/processed/wind.rds")
setup <- readRDS("data/processed/origins_exp2.rds")
ordersCV <- readRDS("results/exp2/orders_cv.rds")
orders1 <- ordersExp1()
h <- setup$horizon
# Smoke test of the pipeline: EXP2_SMOKE=1 uses 2 origins only
if (Sys.getenv("EXP2_SMOKE") == "1") setup$testOrigins <- head(setup$testOrigins, 2)

carimaModel <- function(o, series) {
  force(o); force(series)
  function(s) {
    fit <- fitVariant(s[[series]], o$p, o$d, o$q, o$distribution, o$drift)
    carimaDistribution(fit, series, h)
  }
}
models <- list(persistence = function(s) persistenceDistribution(s, h))
for (m in names(ordersCV)) {
  models[[m]] <- carimaModel(ordersCV[[m]], ordersCV[[m]]$series)
}
for (m in c("auto_carima_logZ", "auto_carima_logZ_cgn", "auto_carima_z", "auto_carima_z_cgn")) {
  models[[m]] <- carimaModel(orders1[[m]], if (grepl("logZ", m)) "logZ" else "z")
}
for (pkg in c("arima", "msarima")) {
  for (pair in c("speed_angle", "logSpeed_angle", "re_im")) {
    models[[paste0(pkg, "_", pair)]] <- local({
      pkgLocal <- pkg; pairLocal <- pair
      function(s) pairDistribution(s, pkgLocal, pairLocal, orders1, h)
    })
  }
}

run <- commandArgs(trailingOnly = TRUE)
if (length(run) == 0) run <- names(models)
stopifnot(all(run %in% names(models)))

empty <- list(point = rep(NA_complex_, h), speed = matrix(NA, h, length(taus)),
              direction = matrix(NA, h, length(taus)))
startTime <- Sys.time()
results <- foreach(origin = setup$testOrigins) %dopar% {
  set.seed(origin)
  s <- windowSeries(wind, (origin - setup$windowLength + 1):origin)
  lapply(setNames(run, run), function(m) tryCatch(models[[m]](s), error = function(e) empty))
}

origins <- setup$testOrigins
for (m in run) {
  out <- list(point = matrix(NA_complex_, length(origins), h, dimnames = list(origins, 1:h)),
              speed = array(NA, c(length(origins), h, length(taus)), dimnames = list(origins, 1:h, taus)),
              direction = array(NA, c(length(origins), h, length(taus)), dimnames = list(origins, 1:h, taus)))
  for (i in seq_along(origins)) {
    out$point[i, ] <- results[[i]][[m]]$point
    out$speed[i, , ] <- results[[i]][[m]]$speed
    out$direction[i, , ] <- results[[i]][[m]]$direction
  }
  saveRDS(out, paste0("results/exp2/forecasts/", m, ".rds"))
  cat(m, ": failed origins", sum(is.na(out$point[, 1])), "\n")
}
actuals <- t(sapply(origins, function(o) wind$z[o + 1:h]))
rownames(actuals) <- origins
saveRDS(actuals, "results/exp2/forecasts/actuals.rds")
cat("models:", length(run), " origins:", length(origins), " time:", format(Sys.time() - startTime), "\n")
