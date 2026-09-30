# Parameters of the selected models, estimated on the selection window (the first eligible
# 4-week window). In the rolling origin they are re-estimated at every origin; these are for
# reporting only.
# cARIMA: each model is fitted from the default start and from nested models (p-1, q) and
# (p, q-1), as in 03, plus the auto.carima estimates for the auto models. logLik is compared with
# carima_order_grid.csv to check that the same optimum was found.
# Output: results/tables/parameters_carima.csv, results/tables/parameters_benchmarks.csv
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
source("R/functions.R")
library(doMC)
registerDoMC(min(parallel::detectCores(), 30))

wind <- readRDS("data/processed/wind.rds")
setup <- readRDS("data/processed/origins.rds")
s <- windowSeries(wind, (setup$selectionOrigin - setup$windowLength + 1):setup$selectionOrigin)
ordersTable <- read.csv("results/tables/orders.csv", stringsAsFactors = FALSE)
orders <- readRDS("results/orders.rds")
grid <- read.csv("results/tables/carima_order_grid.csv", stringsAsFactors = FALSE)

formatComplex <- function(x) sprintf("%.4f%+.4fi", Re(x), Im(x))

# cARIMA(p, d, q) with warm starts from the nested models, memoised by (p, q)
fitNested <- function(y, p, d, q, distribution, drift, extraStarts = list()) {
  cache <- list()
  fitPQ <- function(p, q) {
    k <- paste(p, q)
    if (!is.null(cache[[k]])) return(cache[[k]])
    start <- function(pFrom, qFrom) {
      f <- fitPQ(pFrom, qFrom)
      c(padStart(coef(f), pFrom, qFrom, p, q, drift),
        if (distribution == "dcgnorm") f$other$shape)
    }
    starts <- list()
    if (p > 0) starts$ar <- start(p - 1, q)
    if (q > 0) starts$ma <- start(p, q - 1)
    cache[[k]] <<- fitCARIMA(y, p, q, starts, distribution, drift, d)
    cache[[k]]
  }
  best <- fitPQ(p, q)
  if (length(extraStarts) > 0) {
    other <- fitCARIMA(y, p, q, extraStarts, distribution, drift, d)
    if (logLik(other) > logLik(best)) best <- other
  }
  best
}

carimaRows <- ordersTable[ordersTable$distribution != "", ]
series <- ifelse(grepl("logZ", carimaRows$model), "logZ", "z")
params <- foreach(i = seq_len(nrow(carimaRows)), .combine = rbind) %dopar% {
  r <- carimaRows[i, ]
  o <- orders[[r$model]]
  extra <- if (!is.null(o$coef)) list(c(o$coef, if (r$distribution == "dcgnorm") o$shape)) else list()
  fit <- fitNested(s[[series[i]]], r$p, r$d, r$q, r$distribution, r$drift, extra)
  b <- coef(fit)
  # Residual covariance of (Re, Im), comparable across distributions (for dcgnorm fit$scale holds
  # the scale and pseudo-scale of the generalised normal instead)
  e <- as.vector(residuals(fit))
  e <- e[!is.na(e)]
  S <- cov(cbind(Re(e), Im(e)))
  g <- grid[grid$model == r$model & grid$p == r$p & grid$q == r$q, ]
  data.frame(model = r$model, series = series[i], distribution = r$distribution,
             p = r$p, d = r$d, q = r$q, constant = r$drift,
             intercept = if (r$drift) formatComplex(b[1]) else "",
             ar = paste(formatComplex(b[grepl("^yLag", names(b))]), collapse = "; "),
             ma = paste(formatComplex(b[grepl("^eLag", names(b))]), collapse = "; "),
             shape = if (r$distribution == "dcgnorm") fit$other$shape else NA,
             sdRe = sqrt(S[1, 1]), sdIm = sqrt(S[2, 2]), corReIm = S[1, 2] / sqrt(S[1, 1] * S[2, 2]),
             # Non-circularity: |pseudo-variance| / variance, 0 = circular, 1 = on a line
             nonCircularity = Mod(mean((e - mean(e))^2)) / mean(Mod(e - mean(e))^2),
             logLik = as.numeric(logLik(fit)), AICc = AICc(fit),
             logLikSelection = if (nrow(g) == 1) g$logLik else NA)
}
write.csv(params, "results/tables/parameters_carima.csv", row.names = FALSE)

# Real-valued benchmarks, with the orders (and drift/constant) selected in 03
benchRows <- ordersTable[ordersTable$distribution == "", ]
bench <- foreach(i = seq_len(nrow(benchRows)), .combine = rbind) %dopar% {
  r <- benchRows[i, ]
  y <- s[[sub("^(arima|msarima)_", "", r$model)]]
  b <- if (grepl("^arima_", r$model)) {
    coef(Arima(y, order = c(r$p, 1, r$q), include.drift = r$drift))
  } else {
    coef(msarima(y, orders = list(ar = r$p, i = 1, ma = r$q), lags = 1, constant = r$drift))
  }
  data.frame(model = r$model, p = r$p, d = 1, q = r$q, constant = r$drift,
             coefficients = paste(sprintf("%s=%.4f", names(b), b), collapse = "; "))
}
write.csv(bench, "results/tables/parameters_benchmarks.csv", row.names = FALSE)

print(params[, c("model", "p", "d", "q", "constant", "intercept", "ar", "ma", "shape",
                 "logLik", "logLikSelection")], row.names = FALSE)
print(params[, c("model", "sdRe", "sdIm", "corReIm", "nonCircularity")], row.names = FALSE)
print(bench, row.names = FALSE)
