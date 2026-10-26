# Experiment 2: orders of the grid cARIMA variants chosen by cross-validation (vector RMSE over
# h = 1..12 on the CV origins), and the diagnostic of all orders on the test origins.
# Output: results/exp2/orders_cv.rds, results/exp2/tables/cv_rmse.csv (variant x order),
# results/exp2/tables/allorders_test_<measure>.csv, figures/exp2_allorders_test.png
source("R/exp2_functions.R")

wind <- readRDS("data/processed/wind.rds")
setup <- readRDS("data/processed/origins_exp2.rds")
h <- setup$horizon
actualsFor <- function(origins) t(sapply(origins, function(o) wind$z[o + 1:h]))

measures <- function(file) {
  x <- readRDS(file)
  a <- actualsFor(x$origins)
  errorsFor <- function(fun) sapply(seq_len(dim(x$forecasts)[2]), function(j) {
    f <- x$forecasts[, j, ]
    sqrt(mean(fun(f, a)^2, na.rm = TRUE))
  })
  list(vector = errorsFor(function(f, a) Mod(f - a)),
       speed = errorsFor(function(f, a) Mod(f) - Mod(a)),
       direction = errorsFor(function(f, a) wrapAngle(Arg(f) - Arg(a)) * 180 / pi),
       failed = colSums(is.na(x$forecasts[, , 1])), names = dimnames(x$forecasts)[[2]])
}

# Cross-validation
cvTable <- NULL
ordersCV <- list()
for (variant in gridVariants$model) {
  file <- paste0("results/exp2/allorders/cv_", variant, ".rds")
  if (!file.exists(file)) next
  m <- measures(file)
  cvTable <- rbind(cvTable, data.frame(model = variant, order = m$names, p = ordersGrid$p, q = ordersGrid$q,
                                       vector = m$vector, speed = m$speed, direction = m$direction,
                                       failed = m$failed))
  best <- which.min(m$vector)
  v <- gridVariants[gridVariants$model == variant, ]
  ordersCV[[variant]] <- list(p = ordersGrid$p[best], d = 1, q = ordersGrid$q[best], drift = v$constant,
                              distribution = v$distribution, series = v$series)
}
write.csv(cvTable, "results/exp2/tables/cv_rmse.csv", row.names = FALSE)
saveRDS(ordersCV, "results/exp2/orders_cv.rds")
cat("Orders chosen by CV (vector RMSE, h = 1..12):\n")
print(data.frame(model = names(ordersCV), p = sapply(ordersCV, `[[`, "p"), q = sapply(ordersCV, `[[`, "q"),
                 cvVectorRMSE = sapply(names(ordersCV), function(m) min(cvTable$vector[cvTable$model == m]))),
      row.names = FALSE)

# Diagnostic on the test origins: accuracy of every order
testFiles <- list.files("results/exp2/allorders", pattern = "^test_", full.names = TRUE)
if (length(testFiles) > 0) {
  testTable <- NULL
  for (file in testFiles) {
    m <- measures(file)
    variant <- sub("^test_(.*)\\.rds$", "\\1", basename(file))
    testTable <- rbind(testTable, data.frame(model = variant, order = m$names, p = ordersGrid$p,
                                             q = ordersGrid$q, vector = m$vector, speed = m$speed,
                                             direction = m$direction, failed = m$failed))
  }
  write.csv(testTable, "results/exp2/tables/allorders_test.csv", row.names = FALSE)
  png("figures/exp2_allorders_test.png", width = 1800, height = 600 * length(testFiles), res = 150)
  par(mfrow = c(length(testFiles), 3), mar = c(4.5, 4.5, 3, 1))
  for (variant in unique(testTable$model)) {
    x <- testTable[testTable$model == variant, ]
    for (measure in c("vector", "speed", "direction")) {
      plot(NA, xlim = c(0, 3), ylim = range(x[[measure]]), xlab = "q (MA order)", ylab = paste("RMSE", measure),
           main = paste0(variant, ": ", measure), las = 1, bty = "l", xaxt = "n")
      axis(1, at = 0:3)
      for (p in 0:3) {
        lines(0:3, x[[measure]][x$p == p], type = "b", pch = 15 + p, col = p + 1, lwd = 2)
      }
      legend("topright", legend = paste("p =", 0:3), col = 1:4, pch = 15:18, lwd = 2, bty = "n", cex = 0.8)
    }
  }
  dev.off()
  cat("\nTest origins, vector RMSE (h = 1..12) by order:\n")
  print(xtabs(round(vector, 4) ~ model + order, testTable))
}
