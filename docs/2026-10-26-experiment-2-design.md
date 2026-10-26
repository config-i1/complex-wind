# Experiment 2: design (2026-10-26)

Response to the comments on `docs/2026-09-30-experiment.pdf` (Experiment 1). Status: **design agreed
2026-10-26 (decisions below); scripts not written yet, nothing run**.

## Comments and responses

1. **cAR(0) looks strange: investigate, cross-validate the orders (especially p).** Diagnostic on
   8 four-week windows (the selection window and 7 test windows over the year), cARIMA(p,1,q) for
   log Z without constant, each fit from the default start and from Hannan-Rissanen (keep the best):
   - p > 0 is best on 6 of 8 windows with the normal error and on 4 of 8 with the generalised
     normal. On the selection window itself the normal model now prefers (1,1,1) by 10.7 AICc;
     Experiment 1 chose (1,1,2) there from a sub-optimal fit (logLik 489.9 vs 492.5). So part of
     the cAR(0) is an artefact of the optimiser and of selecting on a single February window.
   - In all (1,1,1) fits the AR and MA coefficients nearly cancel: |phi1 + theta1| = 0.03-0.12.
     (1 - phi B) and (1 + theta B) are almost a common factor, so cARIMA(1,1,1) is close to
     cARIMA(0,1,0). The orders are weakly identified, which is why AICc moves between (0,1,q),
     (p,1,0) and (1,1,1) from window to window.
   - Response: (a) a diagnostic run of all 16 (p,q) on the test origins, reported by p;
     (b) orders selected by rolling-origin cross-validation on separate CV origins;
     (c) all cARIMA fits start from Hannan-Rissanen as well as the default.
2. **Horizon up to 2 hours.** h = 1..12 everywhere (selection, forecasts, evaluation).
3. **Phase plane.** Done: `R/08_phase_plane.R` -> `figures/phase_plane.png`: (1) speed |Z| against
   direction Arg Z, (2) |log Z| against Arg(log Z), with log Z = log speed + i direction (radians in
   (-pi, pi]); all of 2018 plus a two-day path. In the log Z plane the two prevailing directions are two
   tight bands (Arg(log Z) about +25 and -50 degrees); crossings of +-180 degrees in the direction are
   jumps of 2 pi in Im(log Z).
4. **cARIMA in logs with a real intercept.** Already run (run 5, `_re` models: constant
   `y ~ 1 - Im(1)`, i.e. drift in log speed only, none in the direction), after the write-up the
   comments refer to. Normal: angle 38.28 (2nd best) but speed 2.366 and vector 3.644, worse than
   persistence; generalised normal: about the same as without constant. Kept in Experiment 2 with
   h = 12 and intervals.
5. **Prediction intervals.** Quantiles tau = 0.01, 0.02, ..., 0.99 for every model, horizon and
   origin; pinball score for each tau; calibration plots (tau on x, achieved minus tau on y).

## Design

**Data and origins.** As in Experiment 1: 4-week windows, the same 300 test origins (all are
eligible for h = 12, which keeps the results comparable). Parameters re-estimated at each origin.

**Order selection (comment 1).**
- CV origins: 100 origins spread over the eligible ones, none within 12 steps of a test origin.
  Their windows overlap test windows (there is no earlier data), but no CV forecast period
  coincides with a test forecast period.
- For each cARIMA variant, all (p,q) in 0..3 x 0..3 with d = 1 are fitted at every CV origin;
  the selected (p,q) minimises the CV vector RMSE over h = 1..12 (decision).
- Diagnostic: the same 16 orders on the test origins for the two main variants (log Z without
  constant, normal and generalised normal), reported as accuracy against p and q.

**Models.**
| model | series, d | constant | error | orders |
|---|---|---|---|---|
| persistence | - | - | empirical | - |
| carima_logZ_nd, _cgn_nd | log Z, 1 | none | normal, gen. normal | CV |
| carima_logZ_re, _cgn_re | log Z, 1 | real only | normal, gen. normal | CV |
| carima_z, _cgn | Z, 1 | complex | normal, gen. normal | CV |
| auto_carima_z_cgn | Z, d chosen (0) | chosen | gen. normal | auto.carima (as run 4) |
| arima_/msarima_ speed_angle, logSpeed_angle, re_im | two real series, 1 | chosen | normal | auto, d = 1 |

Also kept (decision: keep everything): log Z with a complex constant (`carima_logZ`, `_cgn`, orders
by CV) and `auto_carima_logZ`, `auto_carima_logZ_cgn`, `auto_carima_z` (auto.carima, as run 4).

**Predictive distributions.**
- cARIMA: `forecast.carima()` of complex on branch `intervals` (developed 2026-10-26 for this):
  `interval="prediction"` is analytical: the h-step error is sum psi_j e_{t+h-j}, with conjugate and
  direct variances sigma^2 sum |psi_j|^2 and varsigma^2 sum psi_j^2 (exact complex normal for dcnorm);
  `interval="simulated"` draws `nsim` innovations from the fitted distribution (rcnorm / rcgnorm) and
  combines them with psi (needed for dcgnorm, where sums of generalised normals are not generalised
  normal). Both return the conjugate and direct variances per horizon; "simulated" also the values.
  Checked on 300 simulated cARIMA(1,1,1) series: coverage of 80% / 95% intervals 0.785 / 0.950
  (normal), 0.804 / 0.959 (generalised normal, simulated). Parameter uncertainty is ignored.
- From these to speed and direction:
  - log Z models: speed = exp(Re log Z), direction = Im log Z, both monotone in one part, so their
    quantiles are exp(quantiles of Re) and the quantiles of Im (unwrapped direction, see Measures):
    analytical for dcnorm, from the simulated values for dcgnorm;
  - Z models: speed |Z| and direction Arg Z depend on both parts, so quantiles come from values drawn
    from the joint distribution at each horizon: the complex normal with the returned variances
    (dcnorm), or the simulated values (dcgnorm). nsim = 10000.
- Real pairs: normal predictive distributions from each model (forecast::Arima standard errors;
  smooth::msarima intervals). Speed from the speed model (exp of the quantiles for log speed),
  direction from the angle model; for Re/Im, speed and direction by Monte Carlo from the two
  independent normals at each horizon.
- Persistence: last value plus the empirical h-step changes in the window (speed; wrapped angle),
  so it has a predictive distribution too.

**Measures.** For speed (m/s) and direction (degrees) only (decision), each tau, h, model:
- pinball score averaged over origins; also averaged over tau (an approximation of CRPS / 2);
- calibration: share of actuals below the tau-quantile, minus tau; plots with tau on x;
- point forecasts: RMSE of speed, direction and vector, h = 1..12, as before.
- Direction is circular: each actual angle is taken as the representative (angle + 360 k) closest to
  the forecast median, and the quantiles are on the unwrapped scale around it.

**Cost (32 cores).** CV: 16 orders x 6 variants x 100 origins x ~8 s / 30 workers = about 4.5 h.
Diagnostic on the test origins: 16 x 2 x 300 x 8 s / 30 = about 2.5 h. Final forecasts with
simulation: about 1-1.5 h.

## Decisions (2026-10-26)

- Intervals for cARIMA: analytical and simulated, in the complex package (done, branch `intervals`).
- Interval measures for speed and direction only.
- CV criterion: vector RMSE over h = 1..12.
- Models: keep everything (all Experiment 1 models, with h = 12, CV orders for the grid models).

## Packages (2026-10-26)

complex from branch `intervals` (= `partial-intercept` + prediction intervals + predict fix), installed
from the local clone; smooth from GitHub master 6bd28cd.
