# cARIMA wind forecasting — experiments

Experiments for the paper **"Wind Forecasting with Complex ARIMA"** (Sergey Svetunkov, Pierre Pinson,
Ivan Svetunkov). Idea: treat wind as ONE complex variable (speed + direction) and model it with
complex ARIMA (cARIMA), then compare against real-valued benchmarks (ARIMA per component, VAR/VARIMA).

## Files

- `Wind_Forecasting_with_Complex_ARIMA.pdf` — current paper draft (5 pages; intro + literature notes,
  methodology stub, evaluation = list of candidate datasets). Authors' LaTeX source to be uploaded later;
  until then use interim `paper/wind_carima.tex` + `paper/references.bib` (no LaTeX installed here,
  so it is not compiled locally). When the real source arrives, move interim content into it.
- `Svetunkov---Svetunkov---Complex-Valued-Econometrics.pdf` — monograph (Springer 2024), 165 pp.
  **Chapter 6 "Complex Dynamic Models" = PDF pages 121–144**; Chapter 7 examples (cARIMA on BJsales,
  cARIMAX on Seatbelts) = pp. 145–158. Git-ignored (copyrighted). No `pdftotext` on this box:
  extract with `python3 -c "import pypdf"` (pages 0-indexed).
- `R/` scripts · `data/raw/` downloaded data (git-ignored) · `data/processed/` cleaned series ·
  `results/{fits,forecasts,tables}` · `figures/` · `paper/`.

## Software

R 4.6.1. Packages: `complex` 1.0.3 (`clm()`, `cacf`, `cpacf`, `cscale`/`cdescale`, `clog`, `rcnorm`,
`complex2vec`), `greybox` (`AICc`, `temporaldummy`), `smooth` (`adam`/`msarima` for real ARIMA),
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
- Error term is complex normal with covariance and pseudo-covariance (non-circular) —
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
(shared code in `R/functions.R`).
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
- Not yet: Fourier/diurnal terms, VAR, log-form variants, power, other datasets.

## Working rules

- Scripts in `R/` numbered in run order (`01_data.R`, `02_explore.R`, …); each reads from `data/`
  and writes to `results/`/`figures/`, never modifies raw data.
- Parallel runs: `registerDoMC()` once, outside loops; check nothing else is running first; ≤30
  workers on this 32-core box; pin `OMP_NUM_THREADS=1` etc.
- Paper claims must come from saved results, not from memory of earlier runs.

## Package development: complex generalised normal

Heavy-tailed errors motivate a `dcgnorm` likelihood in `clm()`. Package source:
`/home/config/R/Projects/Packages/complex`, branch `dcgnorm`; plan in `dev/PLAN-dcgnorm.md`,
prototype with checks in `dev/prototype-dcgnorm.R`. Only β is estimated; the scale and pseudo-scale
are concentrated out in closed form. Installed in the user library: `complex` 1.0.3 from branch `dcgnorm`
(2026-09-25; includes the AIC/BIC fix: logLik df = all real parameters, multivariate AICc/BICc).
`dcgnorm` itself is not implemented yet. Remote: git@github.com:openforecast-org/complex.git (no ssh in
this container, so the user does push/fetch).
