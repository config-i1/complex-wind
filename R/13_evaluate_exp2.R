# Experiment 2: evaluation of point forecasts and predictive distributions, h = 1..12.
# Point: RMSE of speed (m/s), direction (degrees) and vector (m/s).
# Distributions (speed and direction): pinball score for each quantile tau = 0.01..0.99 and
# calibration (share of actuals below the tau-quantile, minus tau).
# Direction is circular: each actual direction is aligned to the forecast median (tau = 0.5), i.e.
# taken as median + shortest-arc difference, on the scale of the quantiles.
# Output: results/exp2/tables/{rmse,pinball_by_h,pinball_by_tau,calibration}_*.csv,
# figures/exp2_calibration_{speed,direction}.png, figures/exp2_pinball.png
source("R/exp2_functions.R")

actuals <- readRDS("results/exp2/forecasts/actuals.rds")
files <- setdiff(list.files("results/exp2/forecasts", pattern = "\\.rds$"), "actuals.rds")
models <- sub("\\.rds$", "", files)
order1 <- c("persistence", gridVariants$model, "auto_carima_logZ", "auto_carima_logZ_cgn",
            "auto_carima_z", "auto_carima_z_cgn")
models <- c(intersect(order1, models), sort(setdiff(models, order1)))
fc <- lapply(setNames(models, models), function(m) readRDS(paste0("results/exp2/forecasts/", m, ".rds")))
h <- ncol(actuals)
nTau <- length(taus)
toDegrees <- 180 / pi

#### Point forecasts ####
rmse <- function(fun) {
  t(sapply(fc, function(x) {
    e <- fun(x$point, actuals)
    c(sqrt(colMeans(e^2, na.rm = TRUE)), "h1-12" = sqrt(mean(e^2, na.rm = TRUE)))
  }))
}
pointTables <- list(speed = rmse(function(f, a) Mod(f) - Mod(a)),
                    direction = rmse(function(f, a) wrapAngle(Arg(f) - Arg(a)) * toDegrees),
                    vector = rmse(function(f, a) Mod(f - a)))
for (m in names(pointTables)) {
  colnames(pointTables[[m]])[1:h] <- paste0("h", 1:h)
  write.csv(round(pointTables[[m]], 4), paste0("results/exp2/tables/rmse_", m, ".csv"))
}

#### Predictive distributions ####
# Observed values on the scale of the quantiles: origin x h
observed <- function(x, target) {
  if (target == "speed") return(Mod(actuals))
  median <- x$direction[, , which.min(abs(taus - 0.5))]
  (median + wrapAngle(Arg(actuals) - median)) * toDegrees
}
quantilesOf <- function(x, target) if (target == "speed") x$speed else x$direction * toDegrees

pinball <- list()
calibration <- list()
for (target in c("speed", "direction")) {
  pinball[[target]] <- array(NA, c(length(models), h, nTau), dimnames = list(models, 1:h, taus))
  calibration[[target]] <- array(NA, c(length(models), h, nTau), dimnames = list(models, 1:h, taus))
  for (m in models) {
    y <- observed(fc[[m]], target)
    q <- quantilesOf(fc[[m]], target)
    for (k in seq_len(nTau)) {
      below <- y < q[, , k]
      pinball[[target]][m, , k] <- colMeans((y - q[, , k]) * (taus[k] - below), na.rm = TRUE)
      calibration[[target]][m, , k] <- colMeans(y <= q[, , k], na.rm = TRUE) - taus[k]
    }
  }
  byH <- cbind(apply(pinball[[target]], c(1, 2), mean), "h1-12" = apply(pinball[[target]], 1, mean))
  colnames(byH)[1:h] <- paste0("h", 1:h)
  write.csv(round(byH, 5), paste0("results/exp2/tables/pinball_by_h_", target, ".csv"))
  write.csv(round(apply(pinball[[target]], c(1, 3), mean), 5),
            paste0("results/exp2/tables/pinball_by_tau_", target, ".csv"))
  calibrationTable <- rbind(
    data.frame(horizons = "h1-12", model = models, apply(calibration[[target]], c(1, 3), mean), check.names = FALSE),
    do.call(rbind, lapply(c(1, 6, 12), function(j) {
      data.frame(horizons = paste0("h", j), model = models, calibration[[target]][, j, ], check.names = FALSE)
    })))
  write.csv(calibrationTable, paste0("results/exp2/tables/calibration_", target, ".csv"), row.names = FALSE)
}

#### Plots ####
groups <- list("cARIMA, log Z" = grep("^carima_logZ", models, value = TRUE),
               "cARIMA, Z and auto.carima" = c(grep("^carima_z", models, value = TRUE),
                                               grep("^auto_carima", models, value = TRUE)),
               "Real pairs" = grep("^(arima|msarima)_", models, value = TRUE))
palette8 <- c("#2a6fdb", "#e8590c", "#2f9e44", "#ae3ec9", "#c92a2a", "#1098ad", "#f08c00", "#5c940d")
calibrationPlot <- function(target, file, units) {
  png(file, width = 2100, height = 1500, res = 150)
  par(mfcol = c(2, 3), mar = c(4.5, 4.5, 3, 1))
  for (hSet in list(1:h, c(1))) {
    for (g in names(groups)) {
      ms <- c("persistence", groups[[g]])
      values <- sapply(ms, function(m) colMeans(matrix(calibration[[target]][m, hSet, ], length(hSet))))
      plot(NA, xlim = c(0, 1), ylim = range(c(values, -0.05, 0.05), na.rm = TRUE), las = 1, bty = "l",
           xlab = "nominal quantile", ylab = "achieved minus nominal",
           main = paste0(g, ", ", target, ", ", if (length(hSet) > 1) "h = 1..12" else "h = 1"))
      abline(h = 0, col = "grey60")
      for (j in seq_along(ms)) {
        lines(taus, values[, j], col = if (j == 1) "black" else palette8[(j - 2) %% 8 + 1],
              lwd = if (j == 1) 2.5 else 1.8, lty = if (j == 1) 2 else 1)
      }
      legend("topleft", legend = ms, col = c("black", palette8[(seq_along(ms[-1]) - 1) %% 8 + 1]),
             lwd = 2, lty = c(2, rep(1, length(ms) - 1)), bty = "n", cex = 0.65)
    }
  }
  dev.off()
}
calibrationPlot("speed", "figures/exp2_calibration_speed.png")
calibrationPlot("direction", "figures/exp2_calibration_direction.png")

png("figures/exp2_pinball.png", width = 2100, height = 1300, res = 150)
par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3, 1))
for (target in c("speed", "direction")) {
  for (g in names(groups)) {
    ms <- c("persistence", groups[[g]])
    values <- sapply(ms, function(m) apply(pinball[[target]][m, , , drop = FALSE], 2, mean))
    plot(NA, xlim = c(1, h), ylim = range(values, na.rm = TRUE), las = 1, bty = "l", xlab = "horizon",
         ylab = paste("mean pinball,", if (target == "speed") "m/s" else "degrees"),
         main = paste0(g, ", ", target))
    for (j in seq_along(ms)) {
      lines(1:h, values[, j], col = if (j == 1) "black" else palette8[(j - 2) %% 8 + 1],
            lwd = if (j == 1) 2.5 else 1.8, lty = if (j == 1) 2 else 1)
    }
    legend("topleft", legend = ms, col = c("black", palette8[(seq_along(ms[-1]) - 1) %% 8 + 1]),
           lwd = 2, lty = c(2, rep(1, length(ms) - 1)), bty = "n", cex = 0.65)
  }
}
dev.off()

cat("RMSE h1-12:\n")
print(round(sapply(pointTables, function(x) x[, "h1-12"]), 3))
cat("\nMean pinball over tau, h1-12:\n")
print(round(sapply(pinball, function(x) apply(x, 1, mean)), 4))
cat("\nMean absolute calibration error over tau, h1-12:\n")
print(round(sapply(calibration, function(x) apply(abs(apply(x, c(1, 3), mean)), 1, mean)), 4))
