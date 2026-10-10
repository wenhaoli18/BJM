# Interval-censoring simulations

Reproducible simulations behind the experimental interval-censoring
extension of BJM (`survivalSub(Surv(L, R, type = "interval2") ~ ...)`,
`fitIntervalBJM()`, and the interval-censored branches of `predictRisk()`,
`performancePlot()` and `calibrationPlot()`). This folder is excluded from
the package build (`.Rbuildignore`).

## Running

From the package root:

```bash
Rscript simulations/interval-censoring/run-all.R                 # everything
Rscript simulations/interval-censoring/02-imputation-vs-midpoint.R  # one study
N_REPS=2 Rscript simulations/interval-censoring/run-all.R        # smoke test
```

- `N_REPS` sets the replicates per scenario (default: each script's own,
  listed below); `N_CORES` the parallel workers (default: all but one).
- Replicates run in parallel with `parallel::mclapply()` (forking: macOS or
  Linux). Replicate `r` always uses `set.seed(r)`, so the results do not
  depend on `N_CORES`.
- The package is loaded from the working tree with `pkgload::load_all()`;
  the scripts use internal functions (`imputeEventTime()`, `icphFit()`,
  `landmark_predictions()`, ...).
- A path containing non-ASCII characters needs a UTF-8 locale
  (`LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8`).
- Each script writes `results/<name>.rds` (all replicates, the replicate
  count and the git commit) and `results/<name>.txt` (the printed summary).
  `results/` is not version-controlled.

## Data-generating model (`R/data-generating.R`)

| Part | Truth |
|---|---|
| Covariate | `x ~ N(0, 1)` |
| Event time | proportional hazards, `exp(0.5 x)`; baseline `weibull` H0(t) = 0.05 t^1.5 or `gompertz` H0(t) = 0.08 (e^0.25t − 1) |
| Event type (competing risks) | `D ~ Bernoulli(plogis(−1 + 0.3 T + 0.5 x))`, known once detected |
| Biomarker | `y = 1 + 0.3 t − 0.2 T [+ 0.5 D] + 0.4 x + b0 + b1 t + e`, b0 ~ N(0, 0.5²), b1 ~ N(0, 0.1²), e ~ N(0, 0.3²) |
| Visits | at 0, then every U(gap) years until C ~ U(4, 15); biomarker measured and event checked at each visit |
| Observation | event in (L, R] (last negative visit, detecting visit); right-censored at L if not detected |

Visit gaps of 1–3 years are the "common" design. Gaps of 2–6 years are an
extreme design, where a third to a half of the events fall before the first
visit.

## Studies

| Script | Question | Default reps | Approx. time* |
|---|---|---|---|
| `01-survival-submodel.R` | Does the interval-censored PH model recover the coefficient (with honest SEs) and the baseline, with spline and piecewise baselines? | 100 | < 1 min |
| `02-imputation-vs-midpoint.R` | Bias of f(Y \| T) estimates: midpoint vs `fitIntervalBJM()` (spline, piecewise) vs imputation from the true f(T), paired with the true-T fit | 60 | ~20 min |
| `03-stochastic-em.R` | Burn-in; does imputing with all-true parameters remove bias (is the imputation distribution right)? Stochastic vs Monte Carlo EM (K draws per iteration) | 200 | ~30 min |
| `04-early-hazard.R` | How much does the extrapolated hazard before the first visit matter; spline's own lower tail vs `lower_tail = "weibull"` vs piecewise | 60 | ~25 min |
| `05-competing-risks.R` | Event-type model and f(Y \| T, D) with competing risks: midpoint vs `fitIntervalBJM(form_conditional_cr = )` | 60 | ~5 min |
| `06-performance-measures.R` | Landmark AUC/Brier: IPCW (Yang et al., 2026) and model-based estimates vs the true values; with a misspecified survival sub-model | 40 | ~5 min |

\* Wall-clock estimates with 12 workers on an Apple-silicon MacBook Pro,
extrapolated from a 2-replicate smoke run (`N_REPS=2`, about 12 minutes for
everything); the full suite takes roughly 1.5 hours.

## What was found during development

These numbers were obtained during development with code equivalent to these
scripts (n = 300 per replicate). Rerunning the scripts regenerates them, and
small differences are expected where the designs were unified (for example,
the common follow-up end C ~ U(4, 15)).

- Survival sub-model: coefficients unbiased (0.502 / −0.725 vs 0.5 / −0.7);
  model SEs match the empirical SDs.
- f(Y \| T), T coefficient (truth −0.2), 2–6-year visits, difference from
  the true-T fit: midpoint +0.016 to +0.018. The piecewise baseline left
  +0.005 (its hazard is flat within visit intervals); the spline −0.002
  (n.s.). With 1–3-year visits the spline was within 0.001, the
  piecewise baseline still +0.004 (Weibull, about 4 SE) and the midpoint
  about +0.007.
- The variance of the spline estimates came from the hazard extrapolated
  before the first visit: their error correlated 0.7 with the error of
  H0(1). `lower_tail = "weibull"` cut the variance but biased the Gompertz
  scenario (−0.010), so it is not the default.
- Stochastic EM: the chain settles within about 3 iterations. Imputation
  with all-true parameters is unbiased (−0.0005, t = −0.8). The stochastic
  EM with the true f(T) showed −0.0034 (t = −4) in the extreme Gompertz
  design, and Monte Carlo EM shrank this with K (−0.0021 at K = 5, −0.0015
  at K = 20), so about half of it comes from single-draw iterations.
- Competing risks, 1–3-year visits: all three coefficients unbiased with
  `fitIntervalBJM()` (midpoint: −0.011, +0.006, −0.023). With 2–6-year
  visits there was a residual +0.015 (event-type T coefficient) and +0.012
  (D), against −0.026 and −0.078 for the midpoint.
- Performance, 1–3-year visits: the IPCW Brier score was about 0.03
  against a true 0.13 (only about 5 of 60 cases had a certain status), and
  the IPCW AUC was +0.01 with an RMSE of about 0.05. The model-based
  measures stayed within 0.003 of the truth, also with a survival
  sub-model missing `x`. With 3–9-month visits the IPCW Brier score was
  still about 25% too low.
