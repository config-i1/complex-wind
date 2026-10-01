# Model selection with complex::auto.carima(constant = "auto"), which chooses p, q (0..3),
# d (0 or 1) and the constant among none / real / imaginary / both, on the selection window, for
# Z, log Z and diff(log Z), with the complex normal and generalised normal distributions.
# No forecasting: only the selected models and their log-likelihoods.
# Output: results/tables/select_constants.csv, results/tables/select_constants_ICs_<case>.csv
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
source("R/functions.R")
library(doMC)
registerDoMC(min(parallel::detectCores(), 30))

wind <- readRDS("data/processed/wind.rds")
setup <- readRDS("data/processed/origins.rds")
s <- windowSeries(wind, (setup$selectionOrigin - setup$windowLength + 1):setup$selectionOrigin)
series <- list(z = s$z, logZ = s$logZ, dlogZ = diff(s$logZ))

cases <- expand.grid(series = names(series), distribution = c("dcnorm", "dcgnorm"),
                     stringsAsFactors = FALSE)
cases$case <- paste0(cases$series, ifelse(cases$distribution == "dcgnorm", "_cgn", ""))

# Log-likelihood in the units of Z: the density of Z from that of log Z has the Jacobian
# |d log z / dz|^2 = 1 / |z|^2, i.e. -2 log(speed) per observation. For diff(log Z) the first
# observation is conditioned on, so its term is excluded.
jacobian <- c(z = 0, logZ = -2 * sum(log(Mod(s$z))), dlogZ = -2 * sum(log(Mod(s$z[-1]))))

startTime <- Sys.time()
fits <- foreach(i = seq_len(nrow(cases))) %dopar% {
  t0 <- Sys.time()
  fit <- auto.carima(series[[cases$series[i]]], orders = list(ar = 3, i = 1, ma = 3),
                     constant = "auto", distribution = cases$distribution[i])
  list(fit = fit, time = as.numeric(Sys.time() - t0, units = "mins"))
}

results <- do.call(rbind, lapply(seq_len(nrow(cases)), function(i) {
  fit <- fits[[i]]$fit
  o <- fit$other$orders
  write.csv(fit$ICs, paste0("results/tables/select_constants_ICs_", cases$case[i], ".csv"),
            row.names = FALSE)
  data.frame(case = cases$case[i], series = cases$series[i], distribution = cases$distribution[i],
             p = o[1], d = o[2], q = o[3], constant = fit$constant,
             shape = if (is.null(fit$other$shape)) NA else round(fit$other$shape, 3),
             nobs = nobs(fit), nparam = nparam(fit),
             logLik = round(as.numeric(logLik(fit)), 2),
             logLikZ = round(as.numeric(logLik(fit)) + jacobian[[cases$series[i]]], 2),
             AICcTrimmed = round(min(fit$ICs$AICc), 2), modelsFitted = nrow(fit$ICs),
             minutes = round(fits[[i]]$time, 1))
}))
write.csv(results, "results/tables/select_constants.csv", row.names = FALSE)
print(results, row.names = FALSE)
cat("time:", format(Sys.time() - startTime), "\n")
