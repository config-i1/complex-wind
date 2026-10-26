# Experiment 2: test and cross-validation origins for h = 12.
# Test origins: the 300 of Experiment 1 (all of them are eligible for h = 12).
# CV origins: 100 origins spread evenly over the eligible ones, whose forecast periods do not
# overlap with those of the test origins (|cv - test| >= 12).
# Output: data/processed/origins_exp2.rds
source("R/exp2_functions.R")

wind <- readRDS("data/processed/wind.rds")
setup <- readRDS("data/processed/origins.rds")
h <- horizon2
windowLength <- setup$windowLength

# Eligibility as in R/02_prepare_data.R, with the horizon of 12
countBack <- function(x, n) as.numeric(stats::filter(as.numeric(x), rep(1, n), sides = 1))
okAll <- countBack(!is.na(wind$z), windowLength + h)
okHorizon <- countBack(!is.na(wind$z) & !wind$filledLong, h)
t <- windowLength:(nrow(wind) - h)
eligible <- t[okAll[t + h] == windowLength + h & okHorizon[t + h] == h]
stopifnot(all(setup$origins %in% eligible))

distanceToTest <- sapply(eligible, function(o) min(abs(o - setup$origins)))
candidates <- eligible[distanceToTest >= h & eligible != setup$selectionOrigin]
cvOrigins <- candidates[unique(round(seq(1, length(candidates), length.out = 100)))]

saveRDS(list(testOrigins = setup$origins, cvOrigins = cvOrigins, windowLength = windowLength,
             horizon = h, selectionOrigin = setup$selectionOrigin),
        "data/processed/origins_exp2.rds")
cat("eligible:", length(eligible), " CV candidates:", length(candidates), " CV origins:", length(cvOrigins),
    " from", format(wind$time[cvOrigins[1]]), "to", format(wind$time[tail(cvOrigins, 1)]), "\n")
print(table(format(wind$time[cvOrigins], "%m")))
