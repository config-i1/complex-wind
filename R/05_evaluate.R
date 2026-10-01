# RMSE by horizon for speed, angle and Z, in original units.
# Speed: |Z| (m/s). Angle: Arg(Z) difference wrapped to (-180, 180] degrees.
# Z: modulus of the complex error, |Z - Zhat| (m/s), i.e. the length of the wind vector error.
# Output: results/tables/rmse_{speed,angle,z}.csv. Columns of models without a forecast file are
# kept from the existing tables.
source("R/functions.R")

actuals <- readRDS("results/forecasts/actuals.rds")
files <- setdiff(list.files("results/forecasts", pattern = "\\.rds$"), "actuals.rds")
models <- sub("\\.rds$", "", files)
# Fixed display order: persistence, cARIMA, then the real-valued pairs
displayOrder <- function(models) {
  c(intersect(c("persistence", "carima_logZ", "carima_logZ_cgn", "carima_logZ_nd",
                "carima_logZ_cgn_nd", "carima_logZ_re", "carima_logZ_cgn_re", "carima_z", "carima_z_cgn",
                "auto_carima_logZ", "auto_carima_logZ_cgn", "auto_carima_z", "auto_carima_z_cgn"),
              models), sort(grep("^(arima|msarima)_", models, value = TRUE)))
}
models <- displayOrder(models)
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
  tab <- round(rmseTable(errors[[m]]), 4)
  file <- paste0("results/tables/rmse_", m, ".csv")
  if (file.exists(file)) {
    tabOld <- as.matrix(read.csv(file, row.names = 1, check.names = FALSE))
    tab <- cbind(tabOld[, !colnames(tabOld) %in% models, drop = FALSE], tab)
    tab <- tab[, displayOrder(colnames(tab)), drop = FALSE]
  }
  write.csv(tab, file)
  cat("\nRMSE", m, "\n")
  print(round(t(tab[c(paste0("h", shownH), "h1-36"), ]), 3))
}
missing <- sapply(seq_along(models), function(j) sum(is.na(errors$z[, 1, j])))
names(missing) <- models
cat("\nOrigins with failed forecasts per model:\n"); print(missing)
