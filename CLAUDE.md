# cARIMA wind forecasting — experiments

Experiments for the paper **"Wind Forecasting with Complex ARIMA"** (Sergey Svetunkov, Pierre Pinson,
Ivan Svetunkov). Idea: treat wind as ONE complex variable (speed + direction) and model it with
complex ARIMA (cARIMA), then compare against real-valued benchmarks (ARIMA per component, VAR/VARIMA).

## Files

- `Wind_Forecasting_with_Complex_ARIMA.pdf` — current paper draft (5 pages; intro + literature notes,
  methodology stub, evaluation = list of candidate datasets). LaTeX source not in the folder yet;
  put it in `paper/` when available.
- `Svetunkov---Svetunkov---Complex-Valued-Econometrics.pdf` — monograph (Springer 2024), 165 pp.
  **Chapter 6 "Complex Dynamic Models" = PDF pages 121–144**; Chapter 7 examples (cARIMA on BJsales,
  cARIMAX on Seatbelts) = pp. 145–158. Git-ignored (copyrighted). No `pdftotext` on this box:
  extract with `python3 -c "import pypdf"` (pages 0-indexed).
- `R/` scripts · `data/raw/` downloaded data (git-ignored) · `data/processed/` cleaned series ·
  `results/{fits,forecasts,tables}` · `figures/` · `paper/`.

## Software

R 4.6.1. Packages: `complex` 1.0.2 (`clm()`, `cacf`, `cpacf`, `cscale`/`cdescale`, `clog`, `rcnorm`,
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

## Wind representation (draft Section 3)

- Direction θ in degrees is circular (0 ≡ 360); never model it as a real number.
- Complex wind (draft eq. 2): `Z_t = v_t * exp(1i * θ_t)`, θ in radians, so Re = v·cos θ,
  Im = v·sin θ. That is NOT the zonal/meridional pair unless θ is converted: meteorological θ is the
  direction wind blows FROM, clockwise from north, giving zonal u = -v·sin θ, meridional w = -v·cos θ.
  Pick one convention, state it in the paper, apply it consistently (the model is rotation-equivariant,
  so point forecasts do not depend on it, but interpretation of Re/Im does).
- Draft mentions modelling logs: `log(R e^{iφ}) = log R + iφ` — note φ then wraps; the angle needs
  unwrapping or this form needs care. Not decided yet.
- Complex normal Z ⇒ |Z| Rayleigh, arg Z uniform only in the circular zero-mean case; real wind is
  non-circular with non-zero mean.
- Target variables in the task statement: "power and angle" — i.e. wind power/speed magnitude +
  direction. Clarify with the user whether the modulus is wind speed or turbine power output.

## Candidate data (from draft)

- Kaggle Wind Turbine SCADA (berkerisen/wind-turbine-scada-dataset): 10-min, speed, direction,
  power, theoretical power (Turkey, 2018). Draft's current choice.
- Kaggle "predict-the-wind-speed-at-a-wind-turbine"; Domino reference-project-wind-turbine.
- Overviews: Effenberger & Ludwig (2022); github.com/sltzgs/OpenWindSCADA. ECMWF extraction
  suggested as possibly best.
Record source URL, licence, download date and any cleaning in `data/README.md` when data arrives.

## Experiment design (to agree with the user before coding)

- Models: cARIMA (orders by AICc); benchmarks: independent ARIMA on Re and Im (`smooth::msarima`/
  `forecast::auto.arima`), ARIMA on speed + naive/circular-mean for direction, VAR on (Re, Im)
  (`vars`), random walk / persistence (strong at short horizons for wind), possibly VETS (`legion`).
- Rolling-origin evaluation; horizons relevant to wind (e.g. 10-min data, h = 1…36).
- Errors on the complex scale: vector error |Z - Ẑ| (RMSE/MAE of modulus of error), plus speed
  error and angular error (wrapped to [-180, 180]). Later: probabilistic scores (energy score).
- Seeds fixed; save fits/forecasts as `.rds` in `results/`, tables as CSV.

## Working rules

- Scripts in `R/` numbered in run order (`01_data.R`, `02_explore.R`, …); each reads from `data/`
  and writes to `results/`/`figures/`, never modifies raw data.
- Parallel runs: `registerDoMC()` once, outside loops; check nothing else is running first; ≤30
  workers on this 32-core box; pin `OMP_NUM_THREADS=1` etc.
- Paper claims must come from saved results, not from memory of earlier runs.
