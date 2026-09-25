# RMSE by horizon for speed, angle and Z, in original units.
# Speed: |Z| (m/s). Angle: Arg(Z) difference wrapped to (-180, 180] degrees.
# Z: modulus of the complex error, |Z - Zhat| (m/s), i.e. the length of the wind vector error.
# Output: results/tables/rmse_{speed,angle,z}.csv
source("R/functions.R")

actuals <- readRDS("results/forecasts/actuals.rds")
files <- setdiff(list.files("results/forecasts", pattern = "\\.rds$"), "actuals.rds")
models <- sub("\\.rds$", "", files)
# Fixed display order: persistence, cARIMA, then the real-valued pairs
models <- c(intersect(c("persistence", "carima_logZ", "carima_logZ_cgn", "carima_z", "carima_z_cgn"),
                      models), sort(grep("^(arima|msarima)_", models, value = TRUE)))
forecasts <- lapply(setNames(models, models),
                    function(m) readRDS(paste0("results/forecasts/", m, ".rds")))
h <- ncol(actuals)

# errors[[measure]]: array origin x horizon x model
errorArray <- function(fun) {
  simplify2array(lapply(forecasts, function(f) fun(f, actuals)))
}
errors <- list(
  speed = errorArray(function(f, a) Mod(f) - Mod(a)),
  angle = errorArray(function(f, a) wrapAngle(Arg(f) - Arg(a)) * 180 / pi),
  z = errorArray(function(f, a) Mod(f - a))
)

rmseTable <- function(e) {
  byH <- sqrt(apply(e^2, c(2, 3), mean, na.rm = TRUE))
  dimnames(byH) <- list(paste0("h", 1:h), models)
  rbind(byH, "h1-36" = sqrt(apply(e^2, 3, mean, na.rm = TRUE)))
}

shownH <- c(1, 3, 6, 12, 18, 24, 36)
for (m in names(errors)) {
  tab <- rmseTable(errors[[m]])
  write.csv(round(tab, 4), paste0("results/tables/rmse_", m, ".csv"))
  cat("\nRMSE", m, "\n")
  print(round(t(tab[c(paste0("h", shownH), "h1-36"), ]), 3))
}
missing <- sapply(seq_along(models), function(j) sum(is.na(errors$z[, 1, j])))
names(missing) <- models
cat("\nOrigins with failed forecasts per model:\n"); print(missing)
