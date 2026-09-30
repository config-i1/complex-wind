# Select model orders once, on the first eligible training window (d = 1 for all models, except
# auto.carima).
# cARIMA: p, q in 0..3 by AICc (clm), for each configuration of series, distribution and drift.
# auto.carima: complex::auto.carima chooses d (0 or 1), the constant and p, q (0..3).
# auto.arima and auto.msarima: their own search.
# Output: results/orders.rds, results/tables/orders.csv, results/tables/carima_order_grid.csv,
# results/tables/auto_carima_ICs_<model>.csv
# Usage: Rscript R/03_select_orders.R [model ...]   (no arguments = all models, including the
# real-valued benchmarks, which are selected together under the name "benchmarks").
# Selected models are merged into the existing results; the others are kept as they are.
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
source("R/functions.R")
library(doMC)
nWorkers <- min(parallel::detectCores(), 30)
registerDoMC(nWorkers)

wind <- readRDS("data/processed/wind.rds")
setup <- readRDS("data/processed/origins.rds")
idx <- (setup$selectionOrigin - setup$windowLength + 1):setup$selectionOrigin
s <- windowSeries(wind, idx)

# cARIMA configurations. drift = FALSE means no intercept (formula y ~ -1).
configs <- data.frame(
  model = c("carima_logZ", "carima_logZ_cgn", "carima_z", "carima_z_cgn",
            "carima_logZ_nd", "carima_logZ_cgn_nd"),
  series = c("logZ", "logZ", "z", "z", "logZ", "logZ"),
  distribution = c("dcnorm", "dcgnorm", "dcnorm", "dcgnorm", "dcnorm", "dcgnorm"),
  drift = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
  stringsAsFactors = FALSE)

# auto.carima configurations
autoConfigs <- data.frame(
  model = c("auto_carima_z", "auto_carima_z_cgn", "auto_carima_logZ", "auto_carima_logZ_cgn"),
  series = c("z", "z", "logZ", "logZ"),
  distribution = c("dcnorm", "dcgnorm", "dcnorm", "dcgnorm"),
  stringsAsFactors = FALSE)

run <- commandArgs(trailingOnly = TRUE)
if (length(run) == 0) run <- c(configs$model, autoConfigs$model, "benchmarks")
stopifnot(all(run %in% c(configs$model, autoConfigs$model, "benchmarks")))
configs <- configs[configs$model %in% run, ]
autoConfigs <- autoConfigs[autoConfigs$model %in% run, ]

orders <- if (file.exists("results/orders.rds")) readRDS("results/orders.rds") else list()
gridFile <- "results/tables/carima_order_grid.csv"
gridOld <- if (file.exists(gridFile)) read.csv(gridFile, stringsAsFactors = FALSE) else NULL

# cARIMA grid, fitted by level p + q so that each order can start from its nested models
if (nrow(configs) > 0) {
  grid <- merge(expand.grid(p = 0:3, q = 0:3), configs)
  grid$AICc <- NA
  grid$logLik <- NA
  grid$shape <- NA
  coefs <- list()
  shapes <- list()
  key <- function(model, p, q) paste(model, p, q)
  for (level in 0:6) {
    rows <- which(grid$p + grid$q == level)
    fits <- foreach(i = rows) %dopar% {
      p <- grid$p[i]; q <- grid$q[i]; model <- grid$model[i]
      distribution <- grid$distribution[i]; drift <- grid$drift[i]
      # Warm starts from the nested models; for dcgnorm the nested shape is the last element
      start <- function(pFrom, qFrom) {
        k <- key(model, pFrom, qFrom)
        c(padStart(coefs[[k]], pFrom, qFrom, p, q, drift),
          if (distribution == "dcgnorm") shapes[[k]])
      }
      starts <- list()
      if (p > 0) starts$ar <- start(p - 1, q)
      if (q > 0) starts$ma <- start(p, q - 1)
      fit <- fitCARIMA(s[[grid$series[i]]], p, q, starts, distribution, drift)
      shape <- if (distribution == "dcgnorm") fit$other$shape else 1
      list(AICc = AICc(fit), logLik = as.numeric(logLik(fit)), shape = shape, coef = coef(fit))
    }
    for (j in seq_along(rows)) {
      i <- rows[j]
      grid$AICc[i] <- fits[[j]]$AICc
      grid$logLik[i] <- fits[[j]]$logLik
      grid$shape[i] <- fits[[j]]$shape
      coefs[[key(grid$model[i], grid$p[i], grid$q[i])]] <- fits[[j]]$coef
      shapes[[key(grid$model[i], grid$p[i], grid$q[i])]] <- fits[[j]]$shape
    }
  }

  for (m in configs$model) {
    g <- grid[grid$model == m, ]
    best <- g[which.min(g$AICc), ]
    orders[[m]] <- list(p = best$p, q = best$q, drift = best$drift, distribution = best$distribution,
                        coef = coefs[[key(m, best$p, best$q)]], shape = best$shape)
  }
  grid <- grid[, c("model", "series", "distribution", "drift", "p", "q", "AICc", "logLik", "shape")]
  if (!is.null(gridOld)) {
    if (!"model" %in% names(gridOld)) {
      # Old format (before the drift option): all models had drift
      gridOld$model <- paste0("carima_", gridOld$series, ifelse(gridOld$distribution == "dcgnorm", "_cgn", ""))
      gridOld$drift <- TRUE
    }
    grid <- rbind(gridOld[!gridOld$model %in% configs$model, names(grid)], grid)
  }
  write.csv(grid, gridFile, row.names = FALSE)
}

# auto.carima. Its ICs are on auto.carima's trimmed common sample, so they are kept in separate
# tables and not in carima_order_grid.csv.
if (nrow(autoConfigs) > 0) {
  autoFits <- foreach(i = seq_len(nrow(autoConfigs))) %dopar% {
    selectAutoCarima(s[[autoConfigs$series[i]]], autoConfigs$distribution[i])
  }
  for (i in seq_len(nrow(autoConfigs))) {
    m <- autoConfigs$model[i]
    write.csv(autoFits[[i]]$ICs, paste0("results/tables/auto_carima_ICs_", m, ".csv"))
    orders[[m]] <- autoFits[[i]][names(autoFits[[i]]) != "ICs"]
  }
}

# Real-valued benchmarks
if ("benchmarks" %in% run) {
  for (series in c("speed", "logSpeed", "angle", "re", "im")) {
    orders[[paste0("arima_", series)]] <- selectArima(s[[series]])
    orders[[paste0("msarima_", series)]] <- selectMsarima(s[[series]])
  }
}

saveRDS(orders, "results/orders.rds")
ordersTable <- data.frame(model = names(orders),
                          p = sapply(orders, `[[`, "p"),
                          d = sapply(orders, function(o) if (is.null(o$d)) 1 else o$d),
                          q = sapply(orders, `[[`, "q"),
                          drift = sapply(orders, `[[`, "drift"),
                          distribution = sapply(orders, function(o) if (is.null(o$distribution)) "" else o$distribution),
                          shape = sapply(orders, function(o) if (is.null(o$shape)) NA else round(o$shape, 3)))
# Keep the rows of models that are in the old table but not in results/orders.rds
ordersFile <- "results/tables/orders.csv"
if (file.exists(ordersFile)) {
  ordersOld <- read.csv(ordersFile, stringsAsFactors = FALSE)
  ordersTable <- rbind(ordersOld[!ordersOld$model %in% ordersTable$model, ], ordersTable)
}
write.csv(ordersTable, ordersFile, row.names = FALSE)
print(ordersTable, row.names = FALSE)
