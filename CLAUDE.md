# cARIMA wind forecasting — experiments

Experiments for the paper **"Wind Forecasting with Complex ARIMA"** (Sergey Svetunkov, Pierre Pinson,
Ivan Svetunkov). Idea: treat wind as ONE complex variable (speed + direction) and model it with
complex ARIMA (cARIMA), then compare against real-valued benchmarks (ARIMA per component, VAR/VARIMA).

## Repository

git@github.com:config-i1/complex-wind.git (remote `origin`, branch `main`). There is no ssh in this
container, so the user pushes (`git push -u origin main` the first time).

## Files

- `Wind_Forecasting_with_Complex_ARIMA.pdf` — current paper draft (5 pages; intro + literature notes,
  methodology stub, evaluation = list of candidate datasets). Authors' LaTeX source to be uploaded later;
  until then use interim `paper/wind_carima.tex` + `paper/references.bib`. On the 4-core box TinyTeX
  is installed in `~/.TinyTeX` (2026-09-30): compile with
  `Rscript -e 'tinytex::pdflatex("file.tex")'`, which installs missing LaTeX packages itself. When the real source arrives, move interim content into it.
- `Svetunkov---Svetunkov---Complex-Valued-Econometrics.pdf` — monograph (Springer 2024), 165 pp.
  **Chapter 6 "Complex Dynamic Models" = PDF pages 121–144**; Chapter 7 examples (cARIMA on BJsales,
  cARIMAX on Seatbelts) = pp. 145–158. Git-ignored (copyrighted). No `pdftotext` on this box:
  extract with `python3 -c "import pypdf"` (pages 0-indexed).
- `R/` scripts · `data/raw/` downloaded data (git-ignored) · `data/processed/` cleaned series ·
  `results/{fits,forecasts,tables}` · `figures/` · `paper/` · `docs/` (dated write-ups, e.g.
  `docs/2026-09-30-experiment.tex`/`.pdf` = full description and results of Experiment 1, runs 1–4).

## Software

R 4.6.1 on the 32-core box; R 4.5.2 on the 4-core box (installed 2026-09-30, user library).
Packages (32-core box, 2026-10-26): `complex` 1.0.3 from local branch `intervals` (6e1e118, =
GitHub `partial-intercept` + analytical/simulated prediction intervals for cARIMA, not pushed yet;
`partial-intercept` is needed for the constant types), `smooth` from GitHub master 6bd28cd: `clm()`, `carima()`,
`auto.carima()`, `cacf`, `cpacf`, `cscale`/`cdescale`, `clog`, `rcnorm`, `complex2vec`; `greybox` (`AICc`, `temporaldummy`), `smooth` (`adam`/`msarima` for real ARIMA),
`legion` (vector ETS), `forecast`, `vars` (VAR), `tseries` (ADF/KPSS), `doMC`/`foreach`, `data.table`.
Not installed: `MTS` (VARMA), `scoringRules`, `circular` — install if needed.

## cARIMA essentials (from Ch. 6)

- Model: `(1 - φ1 B - … - φp B^p)(1 - B)^d y_t = β0 + (1 + θ1 B + … + θq B^q) ε_t`, all parameters
  complex. Equivalent to a restricted bivariate VARIMA: 2 real parameters per lag instead of 4
  (the 2×2 coefficient matrix has the rotation-scaling form [[a, -b],[b, a]]). This parsimony vs
  VAR is the core comparison argument of the paper.
- `clm()` implements the "direct" form (6.35) — forecasts `ŷ` in original units, no manual
  un-differencing; missing initial values from differencing are back-extrapolated, so models with
  different `d` are comparable by IC.
- Stationarity/invertibility: roots of the complex AR/MA polynomials outside unit circle.
- Identification: cACF/cPACF (conjugate AND direct, use both) are unreliable for order selection;
  select p, q by AICc/BIC, which requires `loss="likelihood"` (default). cMA needs numeric
  optimisation, so CLS loss is not available with MA terms.
- Differencing: ADF/KPSS on each part, or on the 1-D MDS projection
  (`complex2vec(y) |> dist() |> cmdscale(k=1)`); tests at 1% in the book.
- Scale both parts with `cscale(y, scaling="norm")` before fitting when units differ; undo with
  `cdescale()`. `clm(scaling=...)` also has an internal option (default `"normalisation"`).
- Error term is complex normal with conjugate and direct variances (book terms; the direct one is
  the "pseudo-variance" elsewhere), non-circular when the direct variance is not zero —
  `summary()` prints the 2×2 error covariance matrix. Non-circularity = anisotropy = prevailing wind
  direction; worth reporting.

## clm() usage patterns

```r
library(complex)
fit <- clm(y ~ 1, data, orders = c(p, d, q), subset = 1:nTrain)  # loss="likelihood" default
summary(fit); AICc(fit)
# Pure ARIMA forecasts: newdata only sets the horizon
fc <- predict(fit, newdata = matrix(NA, h, 1))   # fc$mean is complex
plot(fc)
# Order search, as in the book (Ch. 7):
orders <- as.matrix(expand.grid(0:3, 0:1, 0:3)); colnames(orders) <- c("cAR","cI","cMA")
fits <- lapply(seq_len(nrow(orders)), \(i) clm(y ~ 1, data, orders = orders[i, ], subset = 1:nTrain))
best <- fits[[which.min(sapply(fits, AICc))]]
```
Exogenous regressors go in the formula (cARIMAX). Use `summary(fit, stepSize=1e-6)` if Hessian-based
standard errors look odd.

## Wind representation (DECIDED 2026-09-24)

- Modulus = **wind speed** v (m/s), not turbine power. θ = wind direction angle.
- `Z_t = v_t * exp(1i * θ_t)`: Re = v·cos θ, Im = v·sin θ are plain Cartesian coordinates of the
  wind vector, NOT zonal/meridional components. Do not convert to meteorological u/w.
- Angle range **(-π, π]**: degrees to radians, then `θ <- ifelse(θ > pi, θ - 2*pi, θ)`. Agreed jump
  from π to -π is accepted as more convenient than [0, 2π).
- Log form `log Z = log v + iθ` is the main target, with d = 1 and θ unwrapped within each window
  (the ±2π jumps of the wrapped angle would otherwise be huge shocks: 696 in T1.csv). Speed 0 -> 0.1.
- Complex normal Z ⇒ |Z| Rayleigh, arg Z uniform only in the circular zero-mean case; real wind is
  non-circular with non-zero mean.

## Data

Kaggle Wind Turbine SCADA, `data/raw/T1.csv` (git-ignored; details, checksum and gap list in
`data/README.md`). One turbine, 2018, 10-min, 50,530 rows, 2,030 missing slots (three gaps of
more than 3 days). Re-download:
`python -m venv v && v/bin/pip install kagglehub && v/bin/python -c "import kagglehub;print(kagglehub.dataset_download('berkerisen/wind-turbine-scada-dataset'))"`.
Direction is bimodal (about 60 degrees and about 200 degrees). Gap handling (split into segments or impute) is still to be decided.
Other candidates later: ECMWF extraction; Effenberger & Ludwig (2022); github.com/sltzgs/OpenWindSCADA.

## Experiment 1: point forecasts (DECIDED 2026-09-25)

Pipeline: `R/02_prepare_data.R` -> `R/03_select_orders.R` -> `R/04_rolling_origin.R` -> `R/05_evaluate.R`
(shared code in `R/functions.R`). `R/06_parameters.R` re-estimates the selected models on the selection
window for reporting (`results/tables/parameters_{carima,benchmarks}.csv`); it found that
`carima_logZ_nd` (1,1,2) was selected from a sub-optimal fit (logLik 489.91 vs 492.51).
- Data: 10-min T1.csv on a regular grid. Gaps of at most 36 slots filled by linear interpolation of Re/Im
  (training only); the forecast horizon may contain only observed values or fills of at most 3 slots.
  Zero speeds replaced by 0.1 (hack, for logs).
- Target: log Z = log(speed) + i·angle, d = 1 (the model works on diff(log speed) and diff(angle)).
  The angle is UNWRAPPED within each training window (cumulative shortest-arc changes), so its
  differences have no ±2π jumps; forecasts are wrapped back to (-π, π].
- Window 4 weeks (4032 obs), h = 1..36, 300 origins spread evenly over the eligible ones after the
  selection window (Mar–Dec, none in Oct because of gaps). Parameters re-estimated at every origin.
- Orders selected ONCE on the first eligible window, fixed for all origins; d = 1 for all models.
  cARIMA: p, q in 0..3 by AICc. clm's default start can hit a poor local optimum (a nested bigger
  model with lower logLik), so every clm fit also starts from nested-model / selection-window
  estimates and keeps the best logLik (`fitCARIMA()`).
- Models: persistence; cARIMA on log Z and on Z (d = 1), each with dcnorm and dcgnorm
  (`_cgn` suffix; orders selected per distribution; warm starts carry the shape as the last
  element of B); independent ARIMA pairs, each with forecast::auto.arima and
  smooth::auto.msarima: (speed, angle), (log speed, angle), (Re, Im). auto.msarima picks d itself, so
  its p, q (and constant) are selected on diff(y) and refitted with i = 1.
- Measures: RMSE by horizon, original units: speed |Z| (m/s), angle (wrapped difference, degrees),
  Z vector error |Z - Ẑ| (m/s). No MAE, no relative measures, no tests, no predictive distributions yet.
- Forecasts: one file per model in `results/forecasts/<model>.rds` (origin x horizon, complex Z);
  `Rscript R/04_rolling_origin.R model1 model2` recomputes only those models.
  `results/forecasts_run1.rds`, `results/orders_run1.rds` = first run (before dcgnorm), kept for reference.
- Run 2 (2026-09-25): no model clearly beats persistence overall. The drift in the angle equation
  rotates forecasts (median 17° at h36 for normal cARIMA on log Z); dcgnorm cuts it to 5°, and
  cARIMA on log Z with dcgnorm is the best at h1–h3.
- Run 3 (2026-09-25): cARIMA on log Z WITHOUT drift (formula y ~ -1; `_nd` suffix), both
  distributions. Needed package fixes for clm without intercept (MA sign, predict, AR check).
  No drift removes the rotation (median 0.2–0.6° at h36). Normal: best angle RMSE of all
  (38.21 h1–36); dcgnorm: best Z at h1–h6. `03_select_orders.R` / `04_rolling_origin.R` take model
  names as arguments and merge into existing results. `results/orders_run2.rds` = orders before run 3.
- Run 4 (2026-09-30): four `auto_carima_*` models (Z and log Z, dcnorm and dcgnorm). Here
  `complex::auto.carima` (stepwise, `orders = list(ar=3, i=1, ma=3)`, `constant = NULL`) chooses
  **d (0 or 1) and the constant** as well as p, q. This is the one exception to "d = 1 for all models".
  Orders are still selected once on the selection window; its ICs (on auto.carima's trimmed common
  sample) are in `results/tables/auto_carima_ICs_<model>.csv`. `orders()` does not work on carima
  objects; the orders are in `fit$other$orders`. Selected: Z (1,0,2) with constant, both
  distributions; log Z normal (1,1,1) with constant; log Z dcgnorm (0,1,2) without constant (the
  same specification as `carima_logZ_cgn_nd`, near-identical results). Run on the 4-core box
  (2.8 h, 0 failures; persistence recomputed there matches the old tables exactly).
  Results (RMSE h1–36): `auto_carima_z_cgn` has the best Z of all models (3.455 vs 3.567 persistence,
  3.522 arima_re_im), from h12 on, but worse speed (2.352 vs 2.282) — d = 0 reverts to the mean
  wind vector, which shrinks |Ẑ|. `auto_carima_z` (normal): Z 3.509, worst speed (2.466).
  `auto_carima_logZ` (with drift) is poor like `carima_logZ` (angle 40.29, Z 3.987).
- Run 5 (2026-10-01): cARIMA on log Z, d = 1, constant with the REAL part only (drift in log speed,
  none in the direction; formula y ~ 1 - Im(1); `_re` suffix), both distributions, grid search and
  default clm start like runs 1-3. Selected (1,1,1) normal (Re drift 0.00146/step, amplified by
  1/(1-phi1) = 2.4) and (0,1,2) dcgnorm (drift ~0). Normal: angle 38.28 (2nd), but speed 2.366 and
  vector 3.644, worse than persistence and than `_nd`; dcgnorm: about the same as `_cgn_nd`.
  A speed drift does not help. `R/07_select_constants.R` (auto.carima, 4 constant types, no
  forecasts) -> `results/tables/select_constants*.csv`. An HR-start re-run of the grid was tried and
  reverted (results in conversation only); `results/orders_run4.rds` = orders before run 5.
- Not yet: Fourier/diurnal terms, VAR, log-form variants, power, other datasets.

## Experiment 2: h = 12, CV orders, intervals (RUN 2026-10-26)

Design and decisions: `docs/2026-10-26-experiment-2-design.md` (response to the comments on the
09-30 write-up). Phase plane: `R/08_phase_plane.R` -> `figures/phase_plane.png` (levels),
`figures/phase_plane_diff.png` (Delta Z, Delta log Z).
Pipeline `R/09`-`R/13` via `results/exp2/run_orders.sh` then `run_forecasts.sh` (32-core box, 30
workers, ~10 h in total, 0 failed fits). Thread variables MUST be exported in the shell before
`Rscript` (Sys.setenv inside R is too late for OpenBLAS).
- CV orders (vector RMSE h1-12, 100 CV origins): p > 0 in 7 of 8 variants, but all 16 orders lie
  within ~1% (AR/MA near-cancellation); best p = 0 order at most 0.6% worse. Test-origin
  diagnostic (`_nd` variants) agrees: range 0.8%, best p = 0 0.15% worse.
- Test RMSE h1-12, vector: best auto_carima_z_cgn 2.291, carima_z_cgn 2.293, carima_logZ_cgn_nd
  2.295 vs persistence 2.307 (<1% gain). Speed: arima_speed_angle 1.545 best, persistence 1.557.
  Direction: carima_logZ_re 27.105, carima_logZ_nd 27.138 vs persistence 27.385.
- Pinball (mean over tau, h1-12): direction best carima_z_cgn 5.52 vs persistence 5.92 (-7%);
  speed persistence 0.410 ~ arima_speed_angle 0.409, carima_z_cgn 0.412. log Z cARIMA
  distributions are poor (speed 0.465 cgn / 0.50 normal).
- Calibration: model intervals are too WIDE in the centre (S-shaped curves: achieved < nominal
  below 0.5, > above) because the normal/Gaussian-ish spread is driven by heavy tails;
  generalised normal halves it. Best calibrated: persistence bootstrap, carima_z_cgn (direction
  0.022), auto_carima_z_cgn (speed 0.015). log Z cARIMA worst (0.04-0.11).

## Working rules

- Scripts in `R/` numbered in run order (`01_data.R`, `02_explore.R`, …); each reads from `data/`
  and writes to `results/`/`figures/`, never modifies raw data.
- Parallel runs: `registerDoMC()` once, outside loops; check nothing else is running first;
  workers = `min(parallel::detectCores(), 30)`; pin `OMP_NUM_THREADS=1` etc.
- Paper claims must come from saved results, not from memory of earlier runs.

## Package development: complex generalised normal

Heavy-tailed errors motivate a `dcgnorm` likelihood in `clm()`. Package source:
`~/R/Projects/Packages/complex` (branch `dcgnorm` merged into `master`); plan in
`dev/PLAN-dcgnorm.md`, prototype with checks in `dev/prototype-dcgnorm.R`. Only β is estimated; the
conjugate and direct scales are concentrated out in closed form. Since 2026-10-01 (branch
`cgnorm-scale-names`) the `dcgnorm()`/`rcgnorm()` arguments are `scaleConjugate`/`scaleDirect`
(formerly `scale`/`pseudoscale`), and for dcgnorm `fit$scale` is `c(scaleConjugate, scaleDirect)`.
Branch `partial-intercept` (on top of it): partial intercept in `clm()` via the formula, `y~x-Im(1)`
(real part only) / `y~x-Re(1)` (imaginary only); `carima(constant=)` is "none"/"real"/"imaginary"/
"both" (TRUE/FALSE dropped); `auto.carima(constant="auto")` tries all four for every d. In this
project `fitCARIMA(drift=)` takes the type, or TRUE/FALSE for the old `orders.rds` entries
(`constantType()`). Same branch `cgnorm-scale-names`: `clm()` stores argument values in its call, so `summary()`/`vcov()` now work on fits
made inside wrappers such as `fitCARIMA()` (before, they failed with "object 'distribution' not found"). Installed in the user library: `complex` 1.0.3 from branch `dcgnorm`
(2026-09-25; includes the AIC/BIC fix: logLik df = all real parameters, multivariate AICc/BICc).
Branch also has: `dcgnorm`/`rcgnorm` + `clm(distribution=)`, `carima()` (pure cARIMA,
msarima-like interface, `forecast()` method), `auto.carima()` (per-d neighbourhood search, ICs on
common trimmed sample, `search="full"`), stationarity/invertibility penalty in clm's loss, fixes for
clm without intercept. Remote: git@github.com:openforecast-org/complex.git (no ssh in
this container, so the user does push/fetch).
