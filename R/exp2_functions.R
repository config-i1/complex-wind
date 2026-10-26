# Shared code for Experiment 2 (h = 12, cross-validated orders, predictive distributions).
# Design: docs/2026-10-26-experiment-2-design.md
source("R/functions.R")

horizon2 <- 12
taus <- seq(0.01, 0.99, by = 0.01)
nsim <- 10000
dir.create("results/exp2/forecasts", recursive = TRUE, showWarnings = FALSE)
dir.create("results/exp2/allorders", recursive = TRUE, showWarnings = FALSE)
dir.create("results/exp2/tables", recursive = TRUE, showWarnings = FALSE)

# cARIMA variants with d = 1 whose orders are chosen by cross-validation
gridVariants <- data.frame(
  model = c("carima_logZ", "carima_logZ_cgn", "carima_logZ_nd", "carima_logZ_cgn_nd",
            "carima_logZ_re", "carima_logZ_cgn_re", "carima_z", "carima_z_cgn"),
  series = c(rep("logZ", 6), "z", "z"),
  distribution = rep(c("dcnorm", "dcgnorm"), 4),
  constant = c("both", "both", "none", "none", "real", "real", "both", "both"),
  stringsAsFactors = FALSE)
ordersGrid <- expand.grid(p = 0:3, q = 0:3)

# Fit cARIMA from the default and the Hannan-Rissanen starting values (keep the higher likelihood)
fitVariant <- function(y, p, d, q, distribution, constant) {
  fitCARIMA(y, p, q, list(), distribution, constant, d, hr = TRUE)
}

#### Predictive distributions: point forecast of Z and quantiles of speed and direction ####
# Each function returns list(point = complex h, speed = h x 99, direction = h x 99), where the
# direction is in radians on an unwrapped scale (any representative; the evaluation aligns it).

# Quantiles by row of a matrix (h x nsim)
rowQuantiles <- function(x) t(apply(x, 1, quantile, probs = taus, names = FALSE))
# Direction of draws relative to the direction of the point forecast (no jumps at +-pi)
alignedArg <- function(draws, centre) Arg(centre) + wrapAngle(Arg(draws) - Arg(centre))

# cARIMA (clm object): analytical intervals for dcnorm, simulated for dcgnorm (complex >= branch
# intervals: predict.clm(interval = "prediction" / "simulated", nsim, scenarios = TRUE))
carimaDistribution <- function(fit, series, h = horizon2) {
  newdata <- matrix(NA, h, 1)
  if (fit$distribution == "dcnorm") {
    f <- predict(fit, newdata = newdata, interval = "prediction")
    if (series == "logZ") {
      # speed = exp(Re log Z) and direction = Im log Z: quantiles of each part
      z <- qnorm(taus)
      return(list(point = exp(f$mean),
                  speed = exp(outer(Re(f$mean), rep(1, length(taus))) + outer(sqrt(f$variances[, "Re"]), z)),
                  direction = outer(Im(f$mean), rep(1, length(taus))) + outer(sqrt(f$variances[, "Im"]), z)))
    }
    # Z: speed and direction depend on both parts; draws from the complex normal of each horizon
    draws <- t(sapply(1:h, function(j) rcnorm(nsim, f$mean[j], f$scaleConjugate[j], f$scaleDirect[j])))
  } else {
    f <- predict(fit, newdata = newdata, interval = "simulated", nsim = nsim, scenarios = TRUE)
    if (series == "logZ") {
      return(list(point = exp(f$mean), speed = rowQuantiles(exp(Re(f$scenarios))),
                  direction = rowQuantiles(Im(f$scenarios))))
    }
    draws <- f$scenarios
  }
  list(point = f$mean, speed = rowQuantiles(Mod(draws)),
       direction = rowQuantiles(alignedArg(draws, f$mean)))
}

# Real benchmarks: normal predictive distribution (mean, sd) for each horizon
arimaNormal <- function(y, spec, h = horizon2) {
  fit <- Arima(y, order = c(spec$p, 1, spec$q), include.drift = spec$drift)
  # forecast() handles the drift; the interval of +-1 standard deviation gives the standard deviation
  f <- forecast::forecast(fit, h = h, level = 100 * (pnorm(1) - pnorm(-1)))
  list(mean = as.numeric(f$mean), sd = as.numeric(f$upper) - as.numeric(f$mean))
}
msarimaNormal <- function(y, spec, h = horizon2) {
  fit <- msarima(y, orders = list(ar = spec$p, i = 1, ma = spec$q), lags = 1, constant = spec$drift)
  # The interval of +-1 standard deviation (level 68.27%) gives the standard deviation
  f <- forecast(fit, h = h, interval = "prediction", level = pnorm(1) - pnorm(-1))
  list(mean = as.numeric(f$mean), sd = as.numeric(f$upper) - as.numeric(f$mean))
}
normalQuantiles <- function(nd) outer(nd$mean, rep(1, length(taus))) + outer(nd$sd, qnorm(taus))

# Pair of real models: "speed_angle", "logSpeed_angle" or "re_im"
pairDistribution <- function(s, pkg, pair, orders, h = horizon2) {
  fun <- if (pkg == "arima") arimaNormal else msarimaNormal
  spec <- function(series) orders[[paste0(pkg, "_", series)]]
  if (pair == "re_im") {
    re <- fun(s$re, spec("re"), h)
    im <- fun(s$im, spec("im"), h)
    point <- complex(real = re$mean, imaginary = im$mean)
    draws <- t(sapply(1:h, function(j) complex(real = rnorm(nsim, re$mean[j], re$sd[j]),
                                               imaginary = rnorm(nsim, im$mean[j], im$sd[j]))))
    return(list(point = point, speed = rowQuantiles(Mod(draws)),
                direction = rowQuantiles(alignedArg(draws, point))))
  }
  angle <- fun(s$angle, spec("angle"), h)
  if (pair == "speed_angle") {
    speed <- fun(s$speed, spec("speed"), h)
    speedQ <- pmax(normalQuantiles(speed), 0)
    speedPoint <- speed$mean
  } else {
    logSpeed <- fun(s$logSpeed, spec("logSpeed"), h)
    speedQ <- exp(normalQuantiles(logSpeed))
    speedPoint <- exp(logSpeed$mean)
  }
  list(point = complex(modulus = pmax(speedPoint, 0), argument = angle$mean),
       speed = speedQ, direction = normalQuantiles(angle))
}

# Persistence with the bootstrap of its historical errors: the h steps ahead changes of speed and
# of direction (shortest arc) within the window
persistenceDistribution <- function(s, h = horizon2) {
  n <- length(s$speed)
  last <- tail(s$z, 1)
  speedQ <- t(sapply(1:h, function(j) quantile(s$speed[(j + 1):n] - s$speed[1:(n - j)], taus, names = FALSE)))
  directionQ <- t(sapply(1:h, function(j) quantile(wrapAngle(s$angle[(j + 1):n] - s$angle[1:(n - j)]),
                                                   taus, names = FALSE)))
  list(point = rep(last, h), speed = pmax(Mod(last) + speedQ, 0), direction = Arg(last) + directionQ)
}

# Orders of the models selected in Experiment 1 (auto.carima models and real benchmarks)
ordersExp1 <- function() {
  tab <- read.csv("results/tables/orders.csv", stringsAsFactors = FALSE)
  setNames(lapply(seq_len(nrow(tab)), function(i) {
    list(p = tab$p[i], d = tab$d[i], q = tab$q[i],
         drift = if (tab$drift[i] %in% c("TRUE", "FALSE")) as.logical(tab$drift[i]) else tab$drift[i],
         distribution = tab$distribution[i])
  }), tab$model)
}
