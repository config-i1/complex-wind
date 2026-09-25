# Histograms of 10-min changes: diff(log(speed)) and diff(angle), consecutive readings only.
d <- read.csv("data/raw/T1.csv", check.names = FALSE)
names(d) <- c("time", "power", "speed", "theor", "dir")
d$time <- as.POSIXct(d$time, format = "%d %m %Y %H:%M", tz = "UTC")

theta <- d$dir * pi / 180
theta <- ifelse(theta > pi, theta - 2 * pi, theta)   # (-pi, pi]

consecutive <- as.numeric(diff(d$time), units = "mins") == 10
positive <- d$speed[-1] > 0 & d$speed[-nrow(d)] > 0

dLogSpeed <- diff(log(d$speed))[consecutive & positive]
dTheta <- diff(theta)[consecutive]                    # wrapped: contains +-2pi jumps
dThetaTrue <- atan2(sin(dTheta), cos(dTheta))         # shortest-arc change (= unwrapped diff)

logHist <- function(x, main, xlab, breaks) {
  h <- hist(x, breaks = breaks, plot = FALSE)
  k <- h$counts > 0
  plot(h$mids[k], h$counts[k], type = "h", xlim = range(h$breaks), lwd = 3, col = "#2a6fdb", log = "y",
       main = main, xlab = xlab, ylab = "count (log scale)", las = 1, bty = "l")
  mtext(sprintf("n=%d  Q1=%.3f  median=%.3f  Q3=%.3f  sd=%.3f", length(x),
                quantile(x, .25), median(x), quantile(x, .75), sd(x)), cex = 0.7)
}

png("figures/diffs_hist.png", width = 1800, height = 600, res = 150)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
logHist(dLogSpeed, "diff(log(speed))", "change in log speed", 200)
logHist(dTheta, "diff(theta), theta in (-pi, pi]", "radians", seq(-2 * pi, 2 * pi, length.out = 201))
logHist(dThetaTrue, "shortest-arc change (unwrapped diff)", "radians", seq(-pi, pi, length.out = 101))
dev.off()

cat("wrapped |dTheta| > pi:", sum(abs(dTheta) > pi), "\n")
cat("dLogSpeed tails: |x|>1:", sum(abs(dLogSpeed) > 1), " min/max:", round(range(dLogSpeed), 2), "\n")
cat("zero-speed pairs dropped:", sum(consecutive & !positive), "\n")
