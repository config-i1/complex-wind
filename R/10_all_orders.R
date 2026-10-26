# Experiment 2: point forecasts of all cARIMA(p,1,q), p, q in 0..3, for a set of origins.
# Usage: Rscript R/10_all_orders.R <cv|test> [variant ...]   (no variants = all grid variants)
# Output: results/exp2/allorders/<set>_<variant>.rds: list(forecasts = complex array
# origin x order x horizon (forecasts of Z), logLik = origin x order, orders = ordersGrid)
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
source("R/exp2_functions.R")
library(doMC)
registerDoMC(min(parallel::detectCores() - 2, 30))

args <- commandArgs(trailingOnly = TRUE)
set <- match.arg(args[1], c("cv", "test"))
variants <- if (length(args) > 1) args[-1] else gridVariants$model
stopifnot(all(variants %in% gridVariants$model))

wind <- readRDS("data/processed/wind.rds")
setup <- readRDS("data/processed/origins_exp2.rds")
origins <- if (set == "cv") setup$cvOrigins else setup$testOrigins
# Smoke test of the pipeline: EXP2_SMOKE=1 uses 2 origins only
if (Sys.getenv("EXP2_SMOKE") == "1") origins <- head(origins, 2)
h <- setup$horizon

for (variant in variants) {
  v <- gridVariants[gridVariants$model == variant, ]
  startTime <- Sys.time()
  tasks <- expand.grid(order = seq_len(nrow(ordersGrid)), origin = seq_along(origins))
  results <- foreach(i = seq_len(nrow(tasks))) %dopar% {
    origin <- origins[tasks$origin[i]]
    o <- ordersGrid[tasks$order[i], ]
    s <- windowSeries(wind, (origin - setup$windowLength + 1):origin)
    tryCatch({
      fit <- fitVariant(s[[v$series]], o$p, 1, o$q, v$distribution, v$constant)
      fc <- forecastCARIMA(fit, h)
      list(forecast = if (v$series == "logZ") exp(fc) else fc, logLik = as.numeric(logLik(fit)))
    }, error = function(e) list(forecast = rep(NA_complex_, h), logLik = NA))
  }
  forecasts <- array(NA_complex_, c(length(origins), nrow(ordersGrid), h),
                     dimnames = list(origins, paste0("(", ordersGrid$p, ",1,", ordersGrid$q, ")"), 1:h))
  logLiks <- matrix(NA, length(origins), nrow(ordersGrid), dimnames = dimnames(forecasts)[1:2])
  for (i in seq_len(nrow(tasks))) {
    forecasts[tasks$origin[i], tasks$order[i], ] <- results[[i]]$forecast
    logLiks[tasks$origin[i], tasks$order[i]] <- results[[i]]$logLik
  }
  saveRDS(list(forecasts = forecasts, logLik = logLiks, orders = ordersGrid, origins = origins),
          paste0("results/exp2/allorders/", set, "_", variant, ".rds"))
  cat(set, variant, ": failed", sum(is.na(forecasts[, , 1])), "of", length(forecasts[, , 1]),
      " time", format(Sys.time() - startTime), "\n")
}
