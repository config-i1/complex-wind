# Phase plane of the wind: magnitude (x) against angle (y), all 10-min observations of 2018
# 1. Z = v exp(i theta): speed v against direction theta.
# 2. log Z = log v + i theta (theta in (-pi, pi]): modulus |log Z| against argument Arg(log Z).
# 3. Delta Z = Z_t - Z_{t-1}: modulus against argument.
# 4. Delta log Z = Delta log v + i Delta theta (shortest-arc change, as in the models with d = 1):
#    modulus against argument.
# Differences only between consecutive observed slots (no differences across gaps).
# Output: figures/phase_plane.png (levels), figures/phase_plane_diff.png (differences)
wind <- readRDS("data/processed/wind.rds")
wrapAngle <- function(x) x - 2 * pi * ceiling((x - pi) / (2 * pi))
logZ <- complex(real = log(wind$speed), imaginary = wind$theta)
deltaZ <- diff(wind$z)
deltaLogZ <- complex(real = diff(log(wind$speed)), imaginary = wrapAngle(diff(wind$theta)))
level <- !is.na(wind$z)

panel <- function(x, angle, xlab, ylab, main, xlim = range(x)) {
  plot(x, angle, pch = 16, cex = 0.2, col = rgb(0.16, 0.44, 0.86, 0.1),
       xlab = xlab, ylab = ylab, main = main, las = 1, bty = "l",
       xlim = xlim, ylim = c(-180, 180), yaxt = "n")
  axis(2, at = seq(-180, 180, 90), las = 1)
}

png("figures/phase_plane.png", width = 1800, height = 800, res = 150)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1), oma = c(2, 0, 0, 0))
panel(wind$speed[level], wind$theta[level] * 180 / pi, "speed |Z| (m/s)", "direction Arg Z (degrees)", "Z")
panel(Mod(logZ[level]), Arg(logZ[level]) * 180 / pi, "modulus |log Z|", "argument Arg(log Z) (degrees)", "log Z")
mtext("All 10-min observations of 2018. log Z = log speed + i direction (radians, in (-pi, pi]).",
      side = 1, line = 0.5, outer = TRUE, cex = 0.8)
dev.off()

ok <- !is.na(deltaZ)
png("figures/phase_plane_diff.png", width = 1800, height = 800, res = 150)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1), oma = c(2, 0, 0, 0))
panel(Mod(deltaZ[ok]), Arg(deltaZ[ok]) * 180 / pi, "modulus |Delta Z| (m/s)",
      "argument Arg(Delta Z) (degrees)", "Delta Z", xlim = c(0, quantile(Mod(deltaZ[ok]), 0.995)))
panel(Mod(deltaLogZ[ok]), Arg(deltaLogZ[ok]) * 180 / pi, "modulus |Delta log Z|",
      "argument Arg(Delta log Z) (degrees)", "Delta log Z",
      xlim = c(0, quantile(Mod(deltaLogZ[ok]), 0.995)))
mtext(paste0("Differences of consecutive 10-min observations of 2018 (", sum(ok), "). ",
             "Delta log Z = Delta log speed + i Delta direction (shortest arc, radians). ",
             "x axes cut at the 99.5% quantile of the modulus."),
      side = 1, line = 0.5, outer = TRUE, cex = 0.75)
dev.off()
cat("Delta log Z modulus quantiles:\n"); print(quantile(Mod(deltaLogZ[ok]), c(0.5, 0.9, 0.99, 0.995, 1)))
