# BJM 0.3.0

## Breaking changes

* `dynamicPrediction()` has been renamed to `predictRisk()`, to read as a
  verb-first pair with the new `predictLongitudinal()` (below) rather than
  as one specifically-named function alongside a generically-named one.
  Existing code that calls `dynamicPrediction()` positionally or by name
  needs to switch to `predictRisk()`; the argument list, return value, and
  `"predictRisk.BJM"` (formerly `"dynamicPrediction.BJM"`) class are
  otherwise unchanged. `checkBandcountConvergence()`'s `predict_fun`
  argument accepts `predictRisk` in place of the old `dynamicPrediction`.
* The objects returned by `predictLongitudinal()` are now classed
  `"predictLongitudinal.BJM"` (one biomarker) and
  `"predictLongitudinalAll.BJM"` (several biomarkers), formerly
  `"dynamicPredictionBio.BJM"` and `"dynamicPredictionBioAll.BJM"`, so the
  class names match the exported function rather than the internal helper
  behind it. Code that checks `inherits(x, "dynamicPredictionBio.BJM")` (or
  the `All` variant) needs to switch to the new names; the object contents
  and the `print()`/`summary()` output are unchanged.

## New features

* Interval-censored event times -- **experimental**. `survivalSub()` now
  accepts `Surv(L, R, type = "interval2")` outcomes (event known only to lie
  in `(L, R]`, e.g. between two clinic visits) and fits them with a
  proportional hazards model by maximum likelihood (new arguments
  `event_time`, `baseline` and `df`). The baseline is by default a
  Royston-Parmar spline for the log cumulative hazard in log time
  (`baseline = "spline"`, `df = 3`); a piecewise-constant hazard
  (`baseline = "piecewise"`) is also available, but in simulations its
  flat hazard within each visit interval biased the imputed event times,
  and so `fitIntervalBJM()`'s longitudinal fit, when the true hazard was
  increasing.
  Because no event time is known exactly, the longitudinal sub-model
  `f(Y | T)` cannot be fit on the observed events as it is for
  right-censored data; the new `fitIntervalBJM()` fills `T` in by
  stochastic EM and multiple imputation, drawing each subject's `T` from
  `f(Y | T) f(T)` on `(L, R]` -- the same conditional density dynamic
  prediction already evaluates -- and pools the resulting
  `longitudinalSub()` fits with `poolLongitudinalSub()`. Its
  `survival_fit_all` and `long_fit_all` plug into `predictRisk()`,
  `predictLongitudinal()` and `simulateTrajectory()` unchanged. A `strata()`
  term gives each stratum its own baseline hazard, as in a stratified Cox
  model. Competing risks are supported with the event type taken as known
  once the event is detected: `fitIntervalBJM(form_conditional_cr = )`
  draws `T` from `f(Y | T, D) P(D | T) f(T)` and refits the event-type
  model, which uses `T`, on every draw. `performancePlot()` and
  `calibrationPlot()` are not yet supported for interval-censored fits.
  Because an event is only detected at a visit, `predictRisk()` and
  `predictLongitudinal()` condition an interval-censored prediction on
  being event-free at each patient's last visit rather than at
  `prediction_time`, and `predictRisk()` also returns `prob_undetected_1`
  (and, with competing risks, `prob_undetected_2`), the probability that
  the event already happened in between; `simulateTrajectory()` requires
  `prediction_time` to be the last visit. `print()` of an interval-censored fit reports how
  many events fell before the subject's first visit, where the baseline
  hazard is extrapolated rather than estimated (see `?survivalSub`).
  Right-censored fits are unaffected: every prediction
  helper now reads the survival model through internal accessors, and a new
  reference test checks their output is unchanged to 1e-12.

* New `simulateTrajectory()` draws complete futures from the fitted backward
  joint model -- an event time, an event type under competing risks, and
  every biomarker's values at chosen times -- conditional on a patient's
  history up to `prediction_time`, by sampling the event time from its
  posterior, the random effects from theirs, and the biomarkers given both.
  A patient with no biomarker measurements is drawn from their baseline
  covariates alone, as a synthetic patient. Ordinal biomarkers of a
  Gaussian-copula fit are supported. Event times past the last follow-up
  time are reported as event-free through `max_event_time` (by default that
  time) rather than extrapolated. Covariates are carried forward from the
  patient's last row, or follow a path given in `future_covariates` (e.g.
  ascites from year 6 on). The draws agree with `predictRisk()` and
  `predictLongitudinal()` within Monte Carlo error.
* New `plot()` method for `simulateTrajectory()` results draws the
  simulated biomarker trajectories (sampled paths, median and interval
  band, and the history conditioned on; category shares for an ordinal
  biomarker) or the cumulative incidence of each event type from the drawn
  event times.
* New `spaghettiPlot()` draws each subject's observed biomarker trajectory,
  optionally colored by eventual outcome with a smoothed mean per group,
  and with `align = "event"` plots against the time remaining until the
  event (`survival time - time`, so the event is at 0), matching the
  backward model's view of the biomarker.
* New `performancePlot()` evaluates the dynamic risk predictions at a set
  of landmark times: at each landmark it predicts every event-free subject's
  risk from their history so far with `predictRisk()` and plots the
  time-dependent AUC and Brier score of that risk over the prediction
  window, with inverse probability of censoring weighting, separately per
  event type under competing risks. Pass held-out data to estimate
  out-of-sample performance.
* New `calibrationPlot()` checks whether the predicted risks are
  numerically right: at each landmark it groups the event-free subjects by
  predicted risk (deciles by default) and plots each group's mean predicted
  risk against its observed risk in the window (one minus Kaplan--Meier, or
  the Aalen--Johansen cumulative incidence per event type under competing
  risks), with 95% confidence intervals and the diagonal for reference.
* New `plot()` methods for fitted and predicted objects, each returning a
  `ggplot`: `plot(long_fit_all)` draws residual, Q-Q and random-effect
  diagnostics per biomarker (`which = "residuals"`, `"qq"`, `"ranef"`) and
  a heat map of the random-effects correlation across biomarkers
  (`which = "corr"`); `plot(survival_fit_all)` draws a forest plot of the
  hazard ratios (and the competing-risks odds ratios) or, with
  `which = "basehaz"`, the baseline cumulative hazard; and
  `plot(predictLongitudinal(...))` draws each patient's predicted biomarker
  density with the point prediction marked (category probabilities for an
  ordinal biomarker).
* New `cifPlot()` plots the Aalen--Johansen cumulative incidence of each
  event type (one minus Kaplan--Meier for a single event type), with
  pointwise confidence bands and optional stratification by a subject-level
  group. It accepts long-format data and either a numeric censoring code or
  `censor_value = NA` (as in `pbc3$status4`).
* New `imputeLongitudinal()` fills in missing longitudinal biomarker values
  before `longitudinalSub()` runs, instead of relying on complete-case
  analysis. `longitudinalSub()` drops rows missing a biomarker's own
  covariates, and then keeps only subjects who have at least one
  observation of *every* biomarker -- so with irregular/interrupted
  follow-up (missing at random given the observed covariates and
  biomarkers), a subject missing just one of several biomarkers at a visit,
  or missing a biomarker's measurements entirely, is dropped from every
  biomarker's fit, not just the one it is missing. `imputeLongitudinal()`
  fits a MIWAE (Missing data Importance-Weighted AutoEncoder; Mattei &
  Frellsen, 2019) jointly across the supplied biomarkers and draws one or
  more plausible completions via self-normalized importance resampling; the
  completed data can be passed straight into `longitudinalSub()` in place
  of the original `data_fit_all`. A second backend, selected with
  `method = "diffusion"`, fits a conditional denoising diffusion
  probabilistic model (DDPM; Ho, Jain & Abbeel, 2020) instead of a VAE,
  trained with the same observed-entries-only masking principle as the
  MIWAE objective (matching MissDiff, Ouyang et al., 2023), and imputes
  missing cells with a RePaint-style (Lugmayr et al., 2022) reverse-
  diffusion sampler that reproduces observed values exactly and only
  extrapolates the missing ones. Both backends share the same call shape,
  return shape, and `diagnostics` output, so switching between them is a
  one-argument change. Requires the optional `torch` package (listed in
  `Suggests`, not `Imports`, so installing/using BJM without it is
  unaffected); every call into `torch` is fully namespaced. The returned
  `diagnostics` also report a classical single-regression-imputation
  baseline alongside the deep-generative completions, so all of these can
  be compared on a given dataset -- deep generative imputation needs
  enough data to fit reliably, and is not automatically the better choice
  at every sample size.
* New `poolLongitudinalSub()` combines `longitudinalSub()` fits across the
  multiple completed datasets returned by `imputeLongitudinal(..., impute =
  "multiple")`, using Rubin's rules (Rubin, 1987) with the Barnard & Rubin
  (1999) small-sample degrees-of-freedom adjustment. `impute = "single"`
  (the default) fills `data_fit_all` with the across-draw mean of the
  generative model's completions and returns one completed dataset -- fast,
  but fitting `longitudinalSub()` on it treats every imputed value as if it
  had been observed, so the resulting standard errors do not reflect the
  uncertainty from not actually knowing the missing values. Fitting
  `longitudinalSub()` once per completion in `data_fit_all_list` (`impute =
  "multiple"`) and combining the fixed-effect estimates with
  `poolLongitudinalSub()` instead propagates that uncertainty into the
  pooled standard errors, degrees of freedom, and p-values, and reports a
  fraction-of-missing-information (`fmi`) per coefficient so it is visible
  which estimates were affected most by the missingness. `print()`-ing the
  result shows a coefficient table per biomarker alongside the FMI values.
* `longitudinalSub()` now accepts categorical (binary/ordinal) biomarkers
  alongside continuous ones, via a Gaussian-copula extension. Every
  biomarker's type -- `"continuous"` or `"ordinal"` -- is auto-detected from
  whether its response column in `data_fit_all` is a factor, or can be set
  explicitly with the new `biomarker_type` argument (which always takes
  priority over auto-detection when supplied). If every biomarker is
  continuous, fitting is completely unchanged: this is the same
  `nlme::lme()`-based code path as before, byte-for-byte. If at least one
  biomarker is ordinal, it is instead fit as a probit cumulative link mixed
  model via `ordinal::clmm()` (requiring the optional `ordinal` package,
  listed in `Suggests`, not `Imports`), representing each observed category
  as an interval on an underlying continuous latent Gaussian score --
  mathematically the same object as a Gaussian copula linking a discrete
  margin to the model's other (continuous or latent) margins. The shared
  random-effects covariance matrix across every biomarker, continuous and
  ordinal alike, is then re-estimated jointly by an ECM
  (Expectation-Conditional-Maximization) algorithm: one new "inner" E-step
  imputes each ordinal observation's latent score (its truncated-normal
  conditional mean given the fitted thresholds and current random-effect
  prediction) before the existing "outer" E-step/M-step -- unchanged from
  `longitudinalSub()`'s original EM -- re-estimates the covariance matrix.
  This inner-imputation step is a deterministic moment-matching plug-in
  (not a full Bayesian/MCEM draw), a documented v1 simplification.
  `predictRisk()` now also works directly on a fitted mixed-type
  model: its internal conditional-density calculation dispatches to a
  Gaussian-copula-aware variant (`conditionalYTCopula()`/
  `conditionalYDTCopula()`, internal, not exported) whenever
  `long_fit_all$biomarker_type` records at least one ordinal marker. Given
  a patient's observed continuous values and observed ordinal categories,
  this variant factors the joint density as the exact multivariate-normal
  density of the continuous sub-vector times the multivariate-normal box
  probability (`mvtnorm::pmvnorm()`) of the ordinal latent scores falling
  in their category-implied intervals, evaluated at the ordinal block's
  distribution conditional on the continuous observations (standard
  multivariate-normal conditioning). For an all-continuous fit
  (`biomarker_type` is `NULL`), dispatch falls through unchanged to the
  original `conditionalYT()`/`conditionalYDT()`, so no existing
  `predictRisk()` behavior is affected. `predictLongitudinal()`
  (predicting a future biomarker *value*, rather than event risk) now also
  dispatches to copula-aware variants the same way -- both its own
  conditional-density denominator (`conditionalYTCopula()`/
  `conditionalYDTCopula()`, shared with `predictRisk()`) and new
  `conditionalYTBioCopula()`/`conditionalYDTBioCopula()` numerators
  (internal, not exported) -- when predicting either a continuous **or an
  ordinal** biomarker from a mixed-type fit; the other jointly-fit
  biomarkers may freely be continuous or ordinal either way. Predicting the
  future *category* of an ordinal biomarker needs no separate code path:
  `conditionalYTBioCopula()`/`conditionalYDTBioCopula()` already bracket the
  candidate category's row between its cumulative-link thresholds via the
  same `biomarker_type`-driven machinery used for every other ordinal row in
  the joint density, and `predictLongitudinal()` represents an ordinal
  `bio_i`'s `Y_all`/`Y_predict` as integer category codes (`1:K`, in
  threshold order) rather than a numeric grid, with `Y_all` carrying a
  `"category_labels"` attribute giving the matching level names; `Y_density`
  is then that category's predicted probability rather than a density.
  `bandcount3` (the candidate-grid resolution) does not apply to an ordinal
  `bio_i`, whose grid is fixed at its category count, and is ignored in that
  case.
* New `predictLongitudinal()` is now the single entry point for predicting
  future biomarker value(s), replacing the two previously separate
  functions `dynamicPredictionBio()`/`dynamicPredictionBioAll()` (now
  internal, un-exported helpers behind it -- direct calls to either from
  existing scripts will need to switch to `predictLongitudinal()`). Pass a
  `bio_i` naming exactly **one** biomarker to predict just that one, which
  returns a single `predictLongitudinal.BJM` object; pass `bio_i` naming
  **more than one** biomarker, or leave it at the default `NULL` (meaning
  every biomarker in `long_fit_all`), to predict several at once, which
  instead returns a named list of such objects (one per biomarker, named by
  that biomarker's response variable, classed `predictLongitudinalAll.BJM`)
  and computes the bio_i-*independent* pipeline stages (restricting to
  at-risk patients, the survival-side integration grid, and the denominator
  conditional density) only **once**, reusing it across every requested
  biomarker instead of recomputing it once per biomarker -- redundant work
  that is especially costly under the Gaussian-copula path above, where
  that denominator involves `mvtnorm::pmvnorm()` Monte-Carlo evaluations.
  `bio_i` may name continuous and/or ordinal biomarkers, in any
  combination. `bandcount2` (the shared survival grid) is auto-tuned once,
  from a single representative biomarker, when left at its default
  `"auto"`, for the multi-biomarker case; `bandcount3` (each biomarker's
  own candidate-value grid) is then auto-tuned separately per biomarker,
  reusing the already-computed shared stage rather than rebuilding it on
  every doubling round -- except for any ordinal biomarker, whose candidate
  grid is fixed at its category count rather than controlled by
  `bandcount3`, so tuning is skipped for it and the returned `"bandcount3"`
  attribute records `NA` for that biomarker. `checkBandcountConvergence()`'s
  `predict_fun` argument now accepts `predictLongitudinal` (in place of the
  now-internal `dynamicPredictionBio`) for the single-biomarker case.

* `predictLongitudinal()` now builds a continuous biomarker's candidate-value
  grid (`Y_all`) in two passes: a cheap coarse pass (step = the biomarker's
  residual SD) over the old wide range finds where any patient's predictive
  density is not negligible, and the `bandcount3` points are then placed
  there only. The grid used to span a fixed 5 times the observed range on
  either side, so most points (about 90% in a simulated example) fell where
  every density was practically 0. The same `bandcount3` now gives a 2-5
  times finer grid, and the `"auto"` starting value of `bandcount3` was
  lowered from 300 to 100. `Y_all` is still a single grid shared by all
  patients predicted, so the returned structure is unchanged.

## Performance

* `predictRisk()` and `predictLongitudinal()` evaluate each patient's joint
  longitudinal density much faster when there are many biomarkers or many
  visits. The patient's covariance (random effects plus residual error,
  across all biomarkers) is now factored once per patient instead of at
  every survival-time grid point, and when the patient has at least three
  observations per random effect it is factored through the Woodbury
  identity, so the full observations-by-observations matrix is never
  formed or inverted. The predictions are unchanged up to floating-point
  rounding (relative differences around 1e-14 on `pbc3`).

## Bug fixes

* `predictPlot()` and `riskPlot()` set line widths with `linewidth`
  instead of `size`, so they no longer trigger ggplot2's "Using `size`
  aesthetic for lines was deprecated" warning; the plots are unchanged.
  BJM now requires ggplot2 >= 3.4.0, where `linewidth` was introduced.
* Past the last follow-up time, the baseline cumulative hazard used by
  `predictRisk()`/`predictLongitudinal()` now continues from its last
  tabulated value with the slope of the least-squares line through the
  table. It used to switch onto that line itself, which does not pass
  through the last value, so the cumulative hazard jumped at the last
  follow-up time: up in `pbc3` (risks there move by about 1.5%), and down
  whenever the hazard increases over follow-up (e.g. Weibull data), which
  gave the interval spanning the last follow-up time a negative probability
  mass.

* With competing risks, the event-type variable (the response of
  `form_conditional_cr`) may now appear in `long_sub_fixed`, as the
  backward model allows. A patient still at risk necessarily has an unknown
  (`NA`) event type, and that `NA` used to be treated as a missing covariate:
  every such patient was dropped and `predictRisk()`/`predictLongitudinal()`
  failed with "replacement has 1 row, data has 0".

* `predictRisk()`/`predictLongitudinal()` now stop with a clear error when
  no patient in `data_predict_all` is at risk at `prediction_time`, or none
  has a measurement of every biomarker. `predictRisk()` used to return an
  empty result without comment and `predictLongitudinal()` failed with
  "'from' must be a finite number".

* The warning that the integration upper limit was capped (at 20 times the
  last follow-up time) is now only given when some patient's conditional
  survival probability beyond the cap exceeds 1%, and reports how much; it
  used to fire whenever it exceeded 0.01%, e.g. for ordinary simulated
  Weibull data. The cap is also never below `prediction_time + horizon`.

* `riskPlot()` failed to draw ("Discrete values supplied to continuous
  scale") for a patient whose survival time is unknown (`NA`), the usual
  case for a new patient: the vertical line marking the event time is now
  left out then.

* `riskPlot()`/`predictPlot()` drew risks on the biomarker's axis scaled by
  twice the largest biomarker value, which is negative or zero when every
  value is `<= 0` (e.g. a log-scale biomarker such as `serBilir`), flipping
  or collapsing the risk axis. They now scale by twice the largest absolute
  value. Both also failed ("'max' not meaningful for factors") when the
  plotted biomarker was ordinal; its history is now plotted as the category
  codes `1:K` that `predictLongitudinal()` predicts on.

* With an ordinal biomarker, `predictRisk()` and `predictLongitudinal()`
  took each observation's category code from the prediction data's own
  factor levels, not the categories the model was fit with. Dropping unused
  levels (e.g. `droplevels()`), reordering them, or passing a character
  column therefore bracketed the latent score between the wrong thresholds,
  silently: in one check a risk of 0.53 instead of 0.10. Ordinal responses
  are now matched to the fitted categories by label, and a value that is
  not one of them is an error.

* Under competing risks, a survival time entering a biomarker's
  `long_sub_fixed` through a transformation (e.g. `log(years)`,
  `I(years^2)`, `poly(years, 2)`) was evaluated at the first integration
  grid point for the whole integral, so `predictRisk()` and
  `predictLongitudinal()` were wrong for such models (without competing
  risks they were correct). Plain `years` and interactions such as
  `years:year` were not affected.

* `predictLongitudinal()`'s candidate grid for a continuous biomarker spans
  the observed values plus five times their range on either side. With a
  single observation, or all observed values equal, the range is 0, the grid
  was a single point, and that point -- the observed value -- was returned
  as the prediction. The grid now uses at least the biomarker's standard
  deviation in the training data as its span.

* `poolLongitudinalSub()` failed ("no applicable method for 'fixef'") on
  fits with an ordinal biomarker. An ordinal biomarker's thresholds and
  slopes are now pooled with Rubin's rules (with an infinite complete-data
  df, as `ordinal::clmm()`'s inference is asymptotic), and the pooled
  `long_fit_all` carries them where prediction reads them.

* `longitudinalSub(biomarker_type = "ordinal")` failed ("response needs to
  be a factor") when the response was numeric, e.g. a 0/1/2 score. A
  numeric response is now treated as ordered by value.

* `imputeLongitudinal(seed = )` seeded only torch, but the MIWAE backend
  draws its completions with R's random number generator, so results were
  not reproducible. `seed` now seeds both, and restores R's random number
  state afterwards so the caller's own stream is unaffected.

* `imputeLongitudinal()`'s documentation now recommends `impute =
  "multiple"`, and explains that the `"single"` (mean-of-draws) completion
  and row-by-row imputation tend to understate the residual and
  random-effects variances, including the `Sigma_fit` used for prediction.

* The Breslow baseline cumulative hazard was read at the tabulated time
  *nearest* to each integration grid point instead of as the step function
  it is. An event just after `prediction_time` (e.g. at 5.002 for
  `prediction_time = 5`) was then counted before the prediction window
  instead of inside it. With few events per unit of time this matters: in
  `pbc3`, patient 2's one-year risk from year 5 rises from 0.037 to 0.039,
  and the two readings converge to these different values as the grid is
  refined, so the difference is not discretization error.

* `longitudinalSub()` failed with "arguments imply differing number of
  rows" when a variable used only in `long_sub_random` had a missing value.
  Such rows are now dropped, as `lme()` already did.

* `predictRisk()` returned zero risks for a negative `horizon`; it is now
  an error.

* `predictRisk()` and `predictLongitudinal()` now document that
  predictions extrapolate the baseline hazard and the longitudinal
  sub-model beyond the last follow-up time, that ordinal biomarkers add
  small Monte Carlo variation (use `set.seed()` for exact reproducibility),
  and that results depend slightly on which patients are predicted
  together.

* `predictRisk()`, `predictLongitudinal()`, and `dynamicPredictionBio()`
  now drop rows of `data_predict_all` measured after `prediction_time`
  (`time_variable > prediction_time`), with a warning, instead of silently
  conditioning on them. The documentation always required the history to
  stop at `prediction_time`, but nothing enforced it, and the README and
  vignette quick-start examples themselves passed patient 2's full history
  (4 measurements after `year = 5`) to a prediction at `prediction_time = 5`,
  which uses information that would not be available at that time. Those
  examples now subset to `year <= 5`. `predictPlot()` and `riskPlot()`
  already truncated per landmark time and are unaffected.

* `predictRisk()` and `predictLongitudinal()` failed with "wrong sign in
  'by' argument" when a patient's survival time was only slightly after
  `prediction_time` (so the baseline-hazard extrapolation grid ended before
  the last training time). The grid is now only extended when it reaches
  past the last training time.

* `predictRisk()` and `predictLongitudinal()` can now predict for a patient
  whose event time and status are not yet known (`NA`), as for a genuinely
  new patient. Previously a missing status gave "non-conformable arrays"
  and a missing survival time gave "'to' must be a finite number". The
  marginal survival model now uses only the right-hand side of
  `form_marginal_surv`, and a patient with a missing survival time is
  treated as at risk (see also the integration upper limit, below).

* A missing biomarker value (or covariate) in one row of `data_predict_all`
  made that patient's whole prediction `NA`, silently. Each biomarker's
  rows with a missing response or covariate are now dropped before
  predicting, so the prediction conditions on the measurements that were
  actually observed. A patient left with no observation of some biomarker
  still cannot be predicted, and a warning now names them.

* `longitudinalSub()` failed with "no applicable method for 'droplevels'"
  unless the patient id column was a factor (as in `pbc3`); numeric and
  character ids now work too. Relatedly, `predictLongitudinal()` looked
  patients up by a factor id's internal integer codes rather than its
  labels, so with factor ids whose codes differ from their labels (e.g.
  ids `"101"`, `"205"`) it could use another patient's covariance; it now
  matches ids by value.

* Biomarkers can now have different random-effects structures (e.g.
  `list(~ 1 | id, ~ year | id)`), or any random-effects formula, not just
  "all random intercepts" or "all random intercept + slope on
  `time_variable`". Previously any other structure failed at prediction
  time with "non-conformable arguments". The random-effects design is now
  built from each biomarker's own `long_sub_random` formula.

* A competing-risks event-type model that does not include the survival
  time (e.g. `status4 ~ age + sex`) made every prediction `NA`, silently. It
  is now handled as event type not depending on event time. `survivalSub()`
  now also rejects a survival time that enters `form_conditional_cr`
  through a transformation or interaction (e.g. `log(years)`,
  `years:age`), which prediction would otherwise have evaluated
  incorrectly; it must enter as a plain main effect.

* A stratified Cox model (`strata()` in `form_marginal_surv`) failed at
  prediction time. Each patient's marginal survival now uses the baseline
  hazard of their own stratum.

* `predictRisk()` and `predictLongitudinal()` no longer use the predicted
  patient's own survival time. The denominator integrates over every event
  time after `prediction_time`; its grid used to stop at twice the largest
  recorded survival time of the patients being predicted -- their future
  outcome, not available at `prediction_time`. The prediction therefore
  changed with that recorded time (and for a patient whose event came soon
  after `prediction_time` the integral stopped early and the risk was
  inflated several-fold, e.g. 0.70 instead of 0.10). The upper limit is now
  the earliest time by which every at-risk patient's model-based
  probability of still being event-free has fallen below 1e-4, and the
  grid's intervals are spaced to hold equal probability mass (so the wider
  range does not slow convergence in `bandcount2`) and start at
  `prediction_time` rather than half an interval before it. Converged
  predictions are unchanged in the cases checked; at coarse bandcounts
  results move by a few percent. Beyond the last follow-up time the
  baseline hazard is, as before, a linear extrapolation.

* `predictRisk()`'s numerator grid integrated over the prediction window
  plus one extra grid interval (half before `prediction_time`, half after
  `prediction_time + horizon`), inflating risks by roughly
  `1 / bandcount1` -- about 10-18% at `bandcount1 = 10`. It now tiles
  exactly `(prediction_time, prediction_time + horizon]`. This was the main
  reason `bandcount1`/`bandcount2 = "auto"` reported "had not converged" on
  essentially every `pbc3` prediction; auto-tuned risks now converge and
  agree with a very fine grid to about 0.5%.

* `predictLongitudinal()`'s predicted value (the mode of the predictive
  density) snapped to the `bandcount3` grid, so it moved in whole grid steps
  as the grid changed -- 5% off at `bandcount3 = 50` in the package's own
  example -- and `bandcount3 = "auto"` reported non-convergence even on a
  fine grid. The mode is now refined between grid points (a parabola
  through the log-density at the peak), and auto-tuning compares a
  predicted biomarker value on the scale of its predictive standard
  deviation rather than relative to the value itself (which exploded for
  values near 0). The predictive density is also no longer passed through
  the `[0, 1]` clamp meant for risk probabilities, which would have
  flattened any density peak above 1.

* `predictRisk()`'s `risk_prob_1`/`risk_prob_2` and `predictLongitudinal()`'s
  `Y_predict` (and `Y_density`'s columns) are now named by patient id.
  Only patients still at risk at `prediction_time` are predicted, and they
  used to come back as an unnamed vector, so a patient dropped as not at
  risk silently shifted every later value. `print()`/`summary()` show the
  ids. With `horizon <= 0`, `predictRisk()` now returns a zero for every
  at-risk patient (and for both causes under competing risks) instead of a
  single `0` and `risk_prob_2 = NULL`.

* `longitudinalSub()` estimated `Sigma_fit` with a design matrix rebuilt
  from the raw formula on the subjects retained in the joint fit, rather
  than with the basis `lme()` was fit with. For `poly()`, `splines::ns()`,
  `splines::bs()` or `factor()` terms this gave a different basis (or
  factor coding) than the coefficients were estimated on whenever some
  subjects were excluded. The design now reuses `lme()`'s own terms and the
  factor levels of the full fitting data. Consequently such terms in a
  continuous biomarker's `long_sub_fixed` no longer trigger a warning (the
  training basis is reused both here and at prediction time); the warning
  now applies to `long_sub_random` and to ordinal biomarkers' fixed
  effects, where the design is still rebuilt from the data at hand.

* The EM estimate of `Sigma_fit` had no iteration limit (all-continuous
  fits), and the mixed continuous/ordinal ECM stopped silently at 100
  iterations -- before converging on the package's own example, which
  needs about 150. They are now capped at 1000 and 500 iterations, with a
  warning if the cap is reached. A singular random-effects covariance, which
  is replaced by the identity matrix, is now reported with a warning
  instead of a bare "An error occurred" message.

* The baseline cumulative hazard past the last follow-up time is now
  evaluated directly from its linear extrapolation instead of being
  tabulated on a fixed grid of step 0.005, which assumed time was measured
  in years: with time in days, one prediction took about 7 seconds (and
  the table ran to millions of rows); it now takes the same 0.1 seconds as
  in years, and gives the same risk.

* `poolLongitudinalSub()` now uses each coefficient's own complete-data
  degrees of freedom from `lme()` (between-subject covariates such as age
  have far fewer than within-subject time terms) instead of the
  intercept's for every coefficient, which overstated the df of
  between-subject effects; computes the fraction of missing information
  with the pooled Barnard-Rubin df, as in `mice::pool()`, rather than the
  complete-data df; and returns `long_fit_all`, a `longitudinalSub.BJM`
  object built from the pooled point estimates (fixed effects, residual
  variances, and averaged `Sigma_fit`) that can be passed to
  `predictRisk()`/`predictLongitudinal()`.

* `parallel` is no longer imported into the namespace (it is only called
  as `parallel::mclapply()`).

* Predictions no longer depend on the biomarkers' units or break down with
  many observations. The conditional densities behind `predictRisk()` and
  `predictLongitudinal()` carry a factor `det(2 * pi * Sigma)^(-1/2)` that
  shrinks rapidly as the biomarkers' scale or the number of observations
  grows. In the no-competing-risk case a `+ 1e-20` added to the denominator
  to avoid dividing by zero then dominated it: multiplying two biomarkers by
  100 (a change of units) took one `pbc3` patient's risk from 0.35 to
  1e-26. At larger scales `det()` overflowed to `Inf`, giving a risk of 0
  (or `NaN` with competing risks), and a patient with an outlying value
  could underflow every density to 0 the same way. The densities are now
  computed as log densities and exponentiated after a per-patient shift
  that cancels in every ratio, and the `+ 1e-20` is gone: the risk is the
  same in any units (0.3481 from x1 to x1e8 in that example). `pbc3`
  results are unchanged.

* `longitudinalSub()` underestimated the joint random-effects covariance
  `Sigma_fit` whenever some subject was excluded from the joint fit (a
  subject is kept only with at least one complete observation of every
  biomarker -- e.g. a subject never measured for one biomarker is
  dropped). The EM update averaged the retained subjects' `E[b b^T]` but
  divided by the number of subjects in the first biomarker's raw data, and
  EM compounded the shortfall over iterations: with 50 of 169 `pbc3`
  subjects never measured for albumin, the diagonal of `Sigma_fit` came
  out at 7%-62% of its correct value, and predicted risks changed by up to
  a factor of 2 once corrected. Fits in which every subject is retained
  (including all of the package's `pbc3` examples) are unchanged.

* `conditionalYT()` and `conditionalYDT()` -- the internal density functions
  behind `predictRisk()`'s risk-probability output whenever 2+
  longitudinal biomarkers are jointly fit -- computed the joint Gaussian
  density of the stacked observation vector across all markers via a "trace
  trick": reducing the quadratic form `(Y - mu)' Sigma_all_solve (Y - mu)`
  to `sum(diag(...))` of a smaller (number-of-markers-sized) matrix. That
  trace only sums the *diagonal* (same-marker) blocks of the reduced
  matrix, which silently discarded every *cross*-marker contribution of
  `Sigma_all_solve` -- i.e. it implicitly treated the biomarkers as
  conditionally independent given the random effects, even though
  `Sigma_all` is deliberately built from the fitted joint random-effects
  covariance (`Sigma_fit`) specifically to capture correlation *between*
  markers. Since `Sigma_fit` is generally not block-diagonal across
  markers, this under- or over-stated `predictRisk()`'s risk
  probabilities whenever the jointly-fit biomarkers had correlated random
  effects -- the ordinary case for a joint model, not an edge case.
  `dynamicPredictionBio()` was not affected: `conditionalYTBio()`/
  `conditionalYDTBio()` already computed this density directly via
  `mvtnorm::dmvnorm()` on the full covariance. Fixed by computing the
  quadratic form directly on the full stacked mean/observation vectors
  instead of via the trace reduction; verified against an independent
  `mvtnorm::dmvnorm()` reference on a real two-biomarker `pbc3` fit with
  correlated random effects (`test-conditional-density-correctness.R`).
  This changes the numeric value (but not the sign, scale of magnitude, or
  validity) of `predictRisk()`'s output for any existing multi-marker
  fit with correlated random effects; the golden-master regression
  baselines (`testdata/baseline.rds`, `testdata/baseline_noCR.rds`) have
  been regenerated to reflect the corrected values.

* `survivalSub()`'s internal `coxph()` call did not pass `x = TRUE, y =
  TRUE`, so the fitted model did not carry its own design matrix and
  response. Downstream functions that call `survival::basehaz()` or
  `survival::survfit()` on this fit (e.g. `predictRisk()` at
  prediction time) can, in that case, need to reconstruct the model frame
  by re-evaluating the fit's captured call -- but they do so in the
  environment of the fit's *formula*, i.e. the caller's environment, not
  `survivalSub()`'s own execution environment. Since `survivalSub()`
  passes its `data_survival_fitting` argument to `coxph()` by that same
  name, this happened to work whenever the caller's data object was also
  literally named `data_survival_fitting` (as in every example, vignette,
  and test shipped with the package) and failed with `object
  'data_survival_fitting' not found` for any other variable name. Fixed
  by passing `x = TRUE, y = TRUE` to the internal `coxph()` call, so the
  fit is self-contained and this reconstruction is never needed.

* `predictRisk()`/`dynamicPredictionBio()` failed with `"non-conformable
  arrays"` when predicting for a subject who had exactly one longitudinal
  observation to condition on, in a single-biomarker (univariate) model.
  The internal helpers `build_conditional_design()` (`R/conditionalDesign.R`)
  and `process_variance()` (`R/processVariance.R`) both built the residual
  covariance piece of the conditional variance as `diag(Sigma_vector)`,
  where `Sigma_vector` holds one residual variance per observed row for
  that subject. Base R's `diag()` is ambiguous on a length-1 numeric
  vector: instead of returning the intended 1x1 diagonal matrix, it
  interprets the single number as a *dimension* and returns an NxN
  identity matrix, which silently corrupted the dimensions of the
  downstream covariance matrix for exactly this one-observation case (any
  subject with 2+ observations was unaffected). Fixed by calling
  `diag(Sigma_vector, length(Sigma_vector))` in both places, which is
  unambiguous for every length, including 1.

# BJM 0.2.0

## Breaking changes

* All function arguments have been renamed to a single consistent
  `snake_case` convention. The package previously mixed `dot.case`,
  `camelCase`, and `PascalCase` for the same underlying concepts (and
  in a few cases used different names for the same concept in
  different functions). Code that calls these functions using
  positional arguments is unaffected; code that uses named arguments
  must update the names. The renames are:

  | Old name | New name |
  | --- | --- |
  | `data.survival.fitting` | `data_survival_fitting` |
  | `data.fit.all` | `data_fit_all` |
  | `data.plot.all` | `data_plot_all` |
  | `data.predict.all` | `data_predict_all` |
  | `data.predict.all.one` (`predictPlot`) | `data_predict_all_one` |
  | `data.predict.all.pre` (`riskPlot`) | `data_predict_all_pre` |
  | `prediction.time` | `prediction_time` |
  | `formMarginalSurv` | `form_marginal_surv` |
  | `formConditionalCR` | `form_conditional_cr` |
  | `survivalVariableAll` | `survival_variable_all` |
  | `survivalTransFunction` | `survival_trans_function` |
  | `LongSubFixed` | `long_sub_fixed` |
  | `LongSubRandom` | `long_sub_random` |

* Removed the `pbc2` example dataset. It was never referenced by any
  exported function, vignette, or test, and `pbc3` -- the dataset
  actually used throughout the package's examples and tests -- is not
  a duplicate of it (`pbc3` recodes `sex`, log-transforms several
  biomarkers, and adds the competing-risk/transformed-time columns
  `status3`, `status4`, `status5`, and `Tyears1`-`Tyears4`). Code that
  called `data(pbc2)` should switch to `data(pbc3)` and account for
  these differences.
* Removed the exported `print_survivalSub()`, `print_longitudinalSub()`,
  `print_dynamicPrediction()`, and `print_dynamicPredictionBio()`
  functions. Each only forwarded to the corresponding S3 `print.*.BJM`
  method, so calling `print()` on these objects (or letting them print
  at the console) is unaffected; only direct calls to the removed
  `print_*` functions need to switch to `print()`. `print_BJM()` has
  been renamed to `printBJM()`, since it prints a pair of independently
  fit sub-model objects together and does not participate in S3
  dispatch on a single `BJM`-classed object.

## New features

* `survivalSub()`, `longitudinalSub()`, `dynamicPrediction()`, and
  `dynamicPredictionBio()` now return named lists (e.g. `risk_prob_1`/
  `risk_prob_2`, `Y_predict`/`Y_density`/`Y_all`) with the fields
  documented in each function's `@return` block, instead of anonymous
  positional lists. Existing code indexing results with `x[[1]]`,
  `x[[2]]`, etc. continues to work unchanged.
* `survivalSub()`, `longitudinalSub()`, `dynamicPrediction()`, and
  `dynamicPredictionBio()` now validate their arguments up front and
  fail with a specific, actionable message (naming the offending
  argument) instead of a cryptic error from deep inside model-fitting
  or indexing code. `predictPlot()`, `riskPlot()`, and `cmtPlot()` got
  the same treatment.
* New `survivalTrans()` helper builds the `survival_variable_all`/
  `survival_trans_function` pair directly from a vector of cut points
  (e.g. `survivalTrans(c(1, 3, 5, 7))`), instead of requiring two
  hand-written, easy-to-misalign parallel lists.
* Every `data_*_all` argument across the pipeline (`data_fit_all`,
  `data_predict_all`, and now also `data_predict_all_one` in
  `predictPlot()` and `data_predict_all_pre` in `riskPlot()`) accepts a
  single bare `data.frame`, reused for every biomarker, instead of a
  repeated list, when all biomarkers share the same measurement data.
* `longitudinalSub()` now warns when a `long_sub_fixed` formula contains
  `poly()` (in its default orthogonal mode), `splines::ns()`,
  `splines::bs()`, or `factor()`. These terms recompute their
  basis/contrasts from whatever data they are given, but
  `dynamicPrediction()`/`dynamicPredictionBio()` rebuild the design
  matrix from a small, patient-specific slice of data at every point on
  the internal prediction grid, not the data the model was fit on -- so
  the basis silently disagrees with the one used at fitting time
  (wrong predictions) or `model.matrix()` fails outright when there are
  too few distinct values. Use `poly(..., raw = TRUE)`, `I(x^2)`,
  `log()`, `sqrt()`, or other terms that do not depend on the
  surrounding data instead.
* `bandcount1`/`bandcount2`/`bandcount3` (in `dynamicPrediction()`,
  `dynamicPredictionBio()`, `predictPlot()`, and `riskPlot()`) no longer
  need to be chosen by hand: they now default to `"auto"` instead of a
  fixed number. Under `"auto"`, the value is started small and doubled,
  comparing the returned predictions to the previous round, until the
  largest relative change drops below 1%, or 2 doublings have been
  tried (so resolving a bandcount costs at most 3 prediction calls, not
  an open-ended loop). `predictPlot()`/`riskPlot()` resolve their
  `"auto"` bandcount(s) once, using a representative probe call, rather
  than repeating the search on every point in their internal
  `horizon`/landmark-time loop. If a bandcount has still not converged
  after hitting this cap, a warning reports it and the result at the
  largest value tried is returned anyway (not an error). Pass an
  explicit number, as in previous package versions, to skip auto-tuning
  and use a fixed value instead.
* New `checkBandcountConvergence()` helper gives direct, manual control
  over the same doubling check that now runs automatically by default
  (e.g. to use a tighter tolerance, or more doublings, than the
  built-in `"auto"` search): it runs `dynamicPrediction()`/
  `dynamicPredictionBio()` once at the bandcount value(s) you supply and
  once more with those value(s) scaled up (by default, doubled), and
  reports the largest relative change in the returned predictions -- at
  the cost of exactly one extra prediction call.

## Bug fixes

* `riskPlot()` built its internal `data_predict_all` accumulator
  without initializing it first, so it could silently pick up a
  leftover object of the same name from the caller's environment.
  Fixed to initialize it explicitly, matching `predictPlot()`.
* The "wrap a single `data.frame` as a list" convenience documented for
  `longitudinalSub()`, `conditionalYT()`, `conditionalYTBio()`,
  `conditionalYDT()`, `conditionalYDTBio()`, and `process_variance()`
  never actually triggered, because `is.list()` is `TRUE` for
  data frames in R. Fixed the guard in all six places to also check
  `is.data.frame()`.
* `cmtPlot()`'s `id_variable` argument was dead code: three internal
  deduplication steps always looked up a literal column named
  `"id_variable"` instead of the column named by the argument, so a
  custom `id_variable` silently had no effect. Fixed to look up the
  specified column.
* `cmtPlot(condi_time2event = NULL)` crashed with `object 'plot_data'
  not found`, because the fallback that picks the midpoint of
  `time_variable` referenced an undefined variable instead of the
  actual `data_plot_all` argument. Fixed, so `condi_time2event = NULL`
  works as documented.
* `dynamicPrediction()`, `dynamicPredictionBio()`, `predictPlot()`, and
  `riskPlot()` only checked that each element of
  `survival_trans_function` was a function, never that it actually
  worked. A transform that throws an error, or returns a character
  value, a vector of the wrong length, or a non-finite value (e.g.
  `log(x)` at `x <= 0`), previously only surfaced as a cryptic failure
  deep inside the per-patient prediction grid, or, in the non-finite
  case, as silently corrupted predictions with no error at all. Fixed
  by test-calling every transform once, up front, on the supplied
  `prediction_time` and validating its output, before any of the
  (potentially expensive) prediction machinery runs.

## Internal changes

* `longitudinalSub()`'s documentation now includes a worked
  `poly()`/`splines::ns()`/`factor()` example, matched by an equivalent
  example in the package's example scripts, showing that the warning
  described above is a prompt to double-check the fitted basis, not a
  sign that predictions from these terms are wrong.
* Fixed `pbc3`'s documentation, which had copied `@usage`/`@format`
  tags from `pbc2` (`data(pbc2)`, "20 variables") instead of describing
  `pbc3` itself (`data(pbc3)`, 27 variables); the variable-by-variable
  `\describe` list was unaffected and already documented all 27 `pbc3`
  columns correctly.
* Added a `tests/testthat` suite, including golden-master
  characterization tests captured from the package's own `pbc3`
  example data, covering both the competing-risk and
  no-competing-risk code paths.
* Extracted shared per-patient design-matrix construction out of
  `conditionalYT()`/`conditionalYDT()` and
  `conditionalYTBio()`/`conditionalYDTBio()` into internal helpers in
  `R/conditionalDesign.R`.
* Extracted shared at-risk subsetting, integration-grid setup, and
  risk-probability clamping logic out of `dynamicPrediction()` and
  `dynamicPredictionBio()` into internal helpers in
  `R/dynamicPredictionShared.R`.
* No numeric behavior change is intended by the internal changes in
  this release; all changes were verified against golden-master
  baselines captured before refactoring.

# BJM 0.1.0

* Initial CRAN release.
