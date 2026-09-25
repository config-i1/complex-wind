# Prepare the wind series and the forecasting origins.
# Input:  data/raw/T1.csv
# Output: data/processed/wind.rds (regular 10-min grid) and data/processed/origins.rds

windowLength <- 4032    # 4 weeks of 10-min observations
horizon <- 36           # 6 hours
maxGapFill <- 36        # fill runs of at most 36 missing slots (6 hours) by interpolation
maxGapActual <- 3       # only fills of at most 3 slots may be used as actuals in the horizon
speedFloor <- 0.1       # zero speeds replaced by 0.1 (to allow logs)
nOrigins <- 300

d <- read.csv("data/raw/T1.csv", check.names = FALSE)
names(d) <- c("time", "power", "speed", "theor", "dir")
d$time <- as.POSIXct(d$time, format = "%d %m %Y %H:%M", tz = "UTC")

# Regular grid
grid <- data.frame(time = seq(as.POSIXct("2018-01-01 00:00", tz = "UTC"),
                              as.POSIXct("2018-12-31 23:50", tz = "UTC"), by = "10 min"))
wind <- merge(grid, d[, c("time", "speed", "dir")], all.x = TRUE)

# Complex wind in Cartesian form, angle in (-pi, pi]
theta <- wind$dir * pi / 180
theta <- ifelse(theta > pi, theta - 2 * pi, theta)
z <- complex(modulus = wind$speed, argument = theta)

# Fill short gaps by linear interpolation of the real and imaginary parts
isMissing <- is.na(z)
runs <- rle(isMissing)
runEnd <- cumsum(runs$lengths)
runStart <- runEnd - runs$lengths + 1
fillable <- runs$values & runs$lengths <= maxGapFill & runStart > 1 & runEnd < length(z)
filled <- rep(FALSE, length(z))
filledLong <- rep(FALSE, length(z))
for (j in which(fillable)) {
  idx <- runStart[j]:runEnd[j]
  before <- runStart[j] - 1
  after <- runEnd[j] + 1
  w <- (idx - before) / (after - before)
  z[idx] <- (1 - w) * z[before] + w * z[after]
  filled[idx] <- TRUE
  filledLong[idx] <- runs$lengths[j] > maxGapActual
}

# Speed floor, then rebuild z so that modulus and angle are consistent
speed <- Mod(z)
speed[!is.na(speed) & speed < speedFloor] <- speedFloor
theta <- Arg(z)                                  # Arg returns (-pi, pi]
z <- complex(modulus = speed, argument = theta)

wind <- data.frame(time = wind$time, speed = speed, theta = theta, z = z,
                   filled = filled, filledLong = filledLong)

# Eligible origins: no missing values in the training window and the horizon, and no values
# from long interpolations in the horizon. An origin is the index of the last training observation.
countBack <- function(x, n) as.numeric(stats::filter(as.numeric(x), rep(1, n), sides = 1))
okAll <- countBack(!is.na(wind$z), windowLength + horizon)
okHorizon <- countBack(!is.na(wind$z) & !filledLong, horizon)
t <- windowLength:(nrow(wind) - horizon)
eligible <- t[okAll[t + horizon] == windowLength + horizon & okHorizon[t + horizon] == horizon]

# The first eligible window is used for order selection; test origins come after it.
selectionOrigin <- eligible[1]
candidates <- eligible[eligible > selectionOrigin]
# Space the origins evenly over the candidates to cover the whole year and times of day
origins <- candidates[unique(round(seq(1, length(candidates), length.out = nOrigins)))]

saveRDS(wind, "data/processed/wind.rds")
saveRDS(list(selectionOrigin = selectionOrigin, origins = origins,
             windowLength = windowLength, horizon = horizon),
        "data/processed/origins.rds")

cat("grid:", nrow(wind), " missing after fill:", sum(is.na(wind$z)),
    " filled:", sum(filled), " (long:", sum(filledLong), ") floored speeds:", sum(d$speed < speedFloor), "\n")
cat("eligible origins:", length(eligible), " selection window ends:",
    format(wind$time[selectionOrigin]), "\n")
cat("test origins:", length(origins), " from", format(wind$time[origins[1]]),
    "to", format(wind$time[tail(origins, 1)]), "\n")
