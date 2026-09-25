# RMSE by horizon for speed, angle and Z, in original units.
# Speed: |Z| (m/s). Angle: Arg(Z) difference wrapped to (-180, 180] degrees.
# Z: modulus of the complex error, |Z - Zhat| (m/s), i.e. the length of the wind vector error.
# Output: results/tables/rmse_{speed,angle,z}.csv
source("R/functions.R")

forecasts <- readRDS("results/forecasts/forecasts.rds")
forecasts <- forecasts[!sapply(forecasts, inherits, "error")]
models <- colnames(forecasts[[1]]$forecast)
h <- nrow(forecasts[[1]]$forecast)

# errors[[measure]]: array origin x horizon x model
errorArray <- function(fun) {
  sapply(forecasts, function(x) fun(x$forecast, x$actual), simplify = "array") |>
    aperm(c(3, 1, 2))
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
