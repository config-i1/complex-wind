# Select model orders once, on the first eligible training window (d = 1 for all models).
# cARIMA: p, q in 0..3 by AICc (clm), separately for the complex normal and generalised normal. auto.arima and auto.msarima: their own search.
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
grid <- expand.grid(p = 0:3, q = 0:3, series = c("logZ", "z"), distribution = c("dcnorm", "dcgnorm"),
                    stringsAsFactors = FALSE)
grid$AICc <- NA
grid$logLik <- NA
grid$shape <- NA
coefs <- list()
shapes <- list()
key <- function(series, distribution, p, q) paste(series, distribution, p, q)
for (level in 0:6) {
  rows <- which(grid$p + grid$q == level)
  fits <- foreach(i = rows) %dopar% {
    p <- grid$p[i]; q <- grid$q[i]; series <- grid$series[i]; distribution <- grid$distribution[i]
    # Warm starts from the nested models; for dcgnorm the nested shape is the last element
    start <- function(pFrom, qFrom) {
        k <- key(series, distribution, pFrom, qFrom)
        c(padStart(coefs[[k]], pFrom, qFrom, p, q), if (distribution == "dcgnorm") shapes[[k]])
    }
    starts <- list()
    if (p > 0) starts$ar <- start(p - 1, q)
    if (q > 0) starts$ma <- start(p, q - 1)
    fit <- fitCARIMA(s[[series]], p, q, starts, distribution)
    shape <- if (distribution == "dcgnorm") fit$other$shape else 1
    list(AICc = AICc(fit), logLik = as.numeric(logLik(fit)), shape = shape, coef = coef(fit))
  }
  for (j in seq_along(rows)) {
    i <- rows[j]
    grid$AICc[i] <- fits[[j]]$AICc
    grid$logLik[i] <- fits[[j]]$logLik
    grid$shape[i] <- fits[[j]]$shape
    coefs[[key(grid$series[i], grid$distribution[i], grid$p[i], grid$q[i])]] <- fits[[j]]$coef
    shapes[[key(grid$series[i], grid$distribution[i], grid$p[i], grid$q[i])]] <- fits[[j]]$shape
  }
}
write.csv(grid, "results/tables/carima_order_grid.csv", row.names = FALSE)

orders <- list()
for (series in c("logZ", "z")) {
  for (distribution in c("dcnorm", "dcgnorm")) {
    g <- grid[grid$series == series & grid$distribution == distribution, ]
    best <- g[which.min(g$AICc), ]
    name <- paste0("carima_", series, ifelse(distribution == "dcgnorm", "_cgn", ""))
    orders[[name]] <- list(p = best$p, q = best$q, drift = TRUE, distribution = distribution,
                           coef = coefs[[key(series, distribution, best$p, best$q)]],
                           shape = best$shape)
  }
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
