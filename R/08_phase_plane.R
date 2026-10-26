# Phase plane of the wind: magnitude (x) against angle (y)
# 1. Z = v exp(i theta): speed v against direction theta.
# 2. log Z = log v + i theta (theta in (-pi, pi]): modulus |log Z| against argument Arg(log Z).
# The cloud shows all observations of 2018; the line follows the wind over two days.
# Output: figures/phase_plane.png
wind <- readRDS("data/processed/wind.rds")
wind <- wind[!is.na(wind$z), ]
logZ <- complex(real = log(wind$speed), imaginary = wind$theta)

# Two days, ending at the end of the selection window
pathStart <- which(wind$time == as.POSIXct("2018-02-25 14:40", tz = "UTC"))
path <- pathStart + 0:287

panel <- function(x, angle, xlab, ylab, main, xlim = range(x)) {
  plot(x, angle, pch = 16, cex = 0.2, col = rgb(0.16, 0.44, 0.86, 0.1),
       xlab = xlab, ylab = ylab, main = main, las = 1, bty = "l",
       xlim = xlim, ylim = c(-180, 180), yaxt = "n")
  axis(2, at = seq(-180, 180, 90), las = 1)
  # Break the line where the angle jumps between -180 and 180
  pathAngle <- angle[path]
  pathAngle[c(FALSE, abs(diff(pathAngle)) > 180)] <- NA
  lines(x[path], pathAngle, col = "#c0392b", lwd = 1.2)
  points(x[path[1]], angle[path[1]], pch = 19, col = "#c0392b")
}

png("figures/phase_plane.png", width = 1800, height = 800, res = 150)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1), oma = c(2, 0, 0, 0))
panel(wind$speed, wind$theta * 180 / pi, "speed |Z| (m/s)", "direction Arg Z (degrees)", "Z")
panel(Mod(logZ), Arg(logZ) * 180 / pi, "modulus |log Z|", "argument Arg(log Z) (degrees)", "log Z")
mtext(paste0("Blue: all 10-min observations of 2018. Red: ", format(wind$time[path[1]], "%d %b %H:%M"),
             " to ", format(wind$time[tail(path, 1)], "%d %b %H:%M"), " (dot = start). ",
             "log Z = log speed + i direction (radians, in (-pi, pi])."),
      side = 1, line = 0.5, outer = TRUE, cex = 0.8)
dev.off()
