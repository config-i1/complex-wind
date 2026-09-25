# Select model orders once, on the first eligible training window (d = 1 for all models).
# cARIMA: p, q in 0..3 by AICc (clm, complex normal). auto.arima and auto.msarima: their own search.
# Output: results/orders.rds and results/tables/orders.csv
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
source("R/functions.R")
library(doMC)
registerDoMC(16)

wind <- readRDS("data/processed/wind.rds")
setup <- readRDS("data/processed/origins.rds")
idx <- (setup$selectionOrigin - setup$windowLength + 1):setup$selectionOrigin
s <- windowSeries(wind, idx)

# cARIMA grid, fitted by level p + q so that each order can start from its nested models
grid <- expand.grid(p = 0:3, q = 0:3, series = c("logZ", "z"), stringsAsFactors = FALSE)
grid$AICc <- NA
grid$logLik <- NA
coefs <- list()
key <- function(series, p, q) paste(series, p, q)
for (level in 0:6) {
  rows <- which(grid$p + grid$q == level)
  fits <- foreach(i = rows) %dopar% {
    p <- grid$p[i]; q <- grid$q[i]; series <- grid$series[i]
    starts <- list()
    if (p > 0) starts$ar <- padStart(coefs[[key(series, p - 1, q)]], p - 1, q, p, q)
    if (q > 0) starts$ma <- padStart(coefs[[key(series, p, q - 1)]], p, q - 1, p, q)
    fit <- fitCARIMA(s[[series]], p, q, starts)
    list(AICc = AICc(fit), logLik = as.numeric(logLik(fit)), coef = coef(fit))
  }
  for (j in seq_along(rows)) {
    i <- rows[j]
    grid$AICc[i] <- fits[[j]]$AICc
    grid$logLik[i] <- fits[[j]]$logLik
    coefs[[key(grid$series[i], grid$p[i], grid$q[i])]] <- fits[[j]]$coef
  }
}
write.csv(grid, "results/tables/carima_order_grid.csv", row.names = FALSE)

orders <- list()
for (series in c("logZ", "z")) {
  g <- grid[grid$series == series, ]
  best <- g[which.min(g$AICc), ]
  orders[[paste0("carima_", series)]] <- list(p = best$p, q = best$q, drift = TRUE,
                                              coef = coefs[[key(series, best$p, best$q)]])
}

# Real-valued benchmarks
for (series in c("speed", "logSpeed", "angle", "re", "im")) {
  orders[[paste0("arima_", series)]] <- selectArima(s[[series]])
  orders[[paste0("msarima_", series)]] <- selectMsarima(s[[series]])
}

saveRDS(orders, "results/orders.rds")
ordersTable <- data.frame(model = names(orders),
                          p = sapply(orders, `[[`, "p"),
                          d = 1,
                          q = sapply(orders, `[[`, "q"),
                          drift = sapply(orders, `[[`, "drift"))
write.csv(ordersTable, "results/tables/orders.csv", row.names = FALSE)
print(ordersTable, row.names = FALSE)
