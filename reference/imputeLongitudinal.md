# Impute missing longitudinal biomarker values

Fills in missing biomarker values before
[`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
is fit, assuming they are missing at random (MAR): whether a value is
missing depends only on what is observed.

**It is often not needed.**
[`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
fits each biomarker's linear mixed model by maximum likelihood on all of
that biomarker's observed values, which is already valid under MAR; only
the joint covariance of the random effects across biomarkers,
`Sigma_fit`, is estimated from the subjects with some value of every
biomarker. In simulations on `pbc3` (see the "Missing data imputation"
article on the package website), no imputation method made the
predictions closer to those from the complete data than this
complete-case analysis when visits were missing values at random; with
half of the subjects never measured for a biomarker, multiple imputation
with `method = "mixed"` was slightly closer. Imputation can also help
when the missingness depends on variables outside the model, if they are
added to the formulas used here.

Three methods are available, selected with `method`:

- `"mixed"` (default) draws each subject's missing values from the
  multivariate linear mixed model
  [`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
  fits, given everything observed for that subject – all visits and all
  biomarkers – with the model's parameters drawn from an approximation
  to their posterior for each completion (proper multiple imputation;
  see
  [`mixed_model_imputations()`](https://wenhaoli18.github.io/BJM/reference/mixed_model_imputations.md)).
  It needs no extra package, keeps the variation between subjects, and
  was the only method in the simulations that left the random-effects
  variances and the event-time coefficients unbiased.

- `"miwae"` fits a MIWAE (Missing data Importance-Weighted AutoEncoder;
  Mattei & Frellsen, 2019) – a variational autoencoder whose training
  objective only scores reconstruction on the observed entries – jointly
  across all supplied biomarkers, and uses self-normalized importance
  resampling from the fitted decoder to draw one or more plausible
  completions of the missing cells.

- `"diffusion"` fits a conditional denoising diffusion probabilistic
  model (DDPM; Ho, Jain & Abbeel, 2020) over the same per-row biomarker
  vector, trained as in CSDI (Tashiro et al., 2021): in every
  mini-batch, each row's observed entries are split at random into
  targets and conditions, and the network learns to denoise the targets
  given the conditions' values, the covariates and the observation mask,
  so the loss only ever involves observed entries. Missing cells are
  then filled by running the reverse diffusion with the missing entries
  as targets and all observed ones as conditions; observed entries are
  returned unchanged.

The two deep generative methods need the torch package. They complete
each visit from that visit's covariates and observed biomarkers (and,
with `history = TRUE`, a summary of the subject's other visits), and in
the simulations they shrank the estimated variation between subjects,
which every prediction uses; with `history = TRUE` they also biased the
event-time coefficients when subjects lacked a biomarker entirely. The
completed data.frame(s) of any method can be passed directly to
[`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
in place of the original `data_fit_all`.

## Usage

``` r
imputeLongitudinal(
  data_fit_all,
  long_sub_fixed,
  long_sub_random,
  time_variable,
  n_imputations = 5,
  latent_dim = 8,
  hidden_units = c(64, 32),
  epochs = 300,
  importance_samples = 20,
  diffusion_steps = 100,
  method = c("mixed", "miwae", "diffusion"),
  impute = c("single", "multiple"),
  history = FALSE,
  seed = NULL
)
```

## Arguments

- data_fit_all:

  As in
  [`longitudinalSub`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md):
  either a single data.frame reused for every biomarker, or a list with
  one data.frame per biomarker (long format, one row per subject-visit).

- long_sub_fixed:

  As in
  [`longitudinalSub`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md):
  a formula, or list of formulas (one per biomarker), whose left-hand
  side must be a bare column name (no transformation) – the column that
  gets imputed.

- long_sub_random:

  As in
  [`longitudinalSub`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md):
  used only to read off the subject id variable (`~ ... | id`).

- time_variable:

  Name of the visit-time column used, together with the id variable, to
  align rows across biomarkers.

- n_imputations:

  Number of stochastic completions to draw. Also controls how many
  completed datasets are returned when `impute = "multiple"`.

- latent_dim:

  Dimension of the VAE's latent space. Only used when
  `method = "miwae"`.

- hidden_units:

  Integer vector giving the hidden-layer sizes of the network: the
  encoder (and, reversed, the decoder) when `method = "miwae"`, or the
  denoising network's body when `method = "diffusion"`.

- epochs:

  Number of training epochs: full-batch gradient steps for
  `method = "miwae"`, passes over the data in mini-batches of 128 rows
  for `method = "diffusion"`. The default suits data sets of a few
  hundred to a few thousand rows; many more epochs can overfit, which
  shows as imputations that vary too little between draws (and so pooled
  standard errors that are too small) – in a simulation, MIWAE's spread
  between draws fell from about the true residual SD at 300 epochs to
  80\\ of it at 3000.

- importance_samples:

  Number of importance samples (`K`) used both in the training objective
  and at imputation time. Only used when `method = "miwae"`.

- diffusion_steps:

  Number of forward/reverse diffusion timesteps. Only used when
  `method = "diffusion"`; more steps generally improve sample quality at
  proportionally higher training/sampling cost.

- method:

  `"mixed"` (default), `"miwae"` or `"diffusion"`; see Description.
  `latent_dim`, `hidden_units`, `epochs`, `importance_samples` and
  `diffusion_steps` only apply to the deep generative methods.

- impute:

  Either `"single"` (default; returns one completed `data_fit_all`,
  filling each missing cell with the mean of the `n_imputations` draws –
  a drop-in replacement for the original `data_fit_all`) or `"multiple"`
  (also returns `n_imputations` separately-drawn completed datasets in
  `data_fit_all_list`, to fit
  [`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
  once per completion and pool the fits with
  [`poolLongitudinalSub`](https://wenhaoli18.github.io/BJM/reference/poolLongitudinalSub.md)).
  Prefer `"multiple"` unless only a little data is missing. The mean of
  several draws is less variable than the values it stands in for, so a
  model fit to the `"single"` completion treats the imputed cells as
  exactly known: its standard errors are too small, and its residual and
  random-effects variances – including `Sigma_fit`, which
  [`predictRisk`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
  and
  [`predictLongitudinal`](https://wenhaoli18.github.io/BJM/reference/predictLongitudinal.md)
  use directly – tend to be underestimated, more so the larger the share
  of imputed cells.

- history:

  For the deep generative methods: whether each visit's inputs include a
  summary of the same subject's other visits (see
  [`history_features()`](https://wenhaoli18.github.io/BJM/reference/history_features.md)).
  Without it an imputed value does not carry the subject's own level,
  which shrinks the estimated variation between subjects; with it, that
  shrinkage is smaller when visits miss values at random, but in
  simulations it biased the event-time coefficients when subjects lacked
  a biomarker entirely. Defaults to `FALSE`.

- seed:

  Optional integer seed for reproducibility. It seeds both torch and R's
  random number generator for the duration of the call; R's random
  number state from before the call is restored afterwards.

## Value

A named list with elements:

- data_fit_all:

  A list with one completed data frame per biomarker (also when
  `data_fit_all` was a single data frame), with every missing biomarker
  cell filled. Each holds one row per (id, `time_variable`) pair seen
  for any biomarker, and only the columns its model uses: the id,
  `time_variable`, the biomarker and the covariates of its
  `long_sub_fixed` formula – other columns of the input are not carried
  over. Ready to pass straight to
  [`longitudinalSub`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
  with the same formulas.

- data_fit_all_list:

  When `impute = "multiple"`, a list of `n_imputations` completed
  versions of `data_fit_all`, each from its own draw; `NULL` when
  `impute = "single"`.

- diagnostics:

  A named list, one element per biomarker (its `deep_generative` entry
  describes the completions of whichever `method` was used), each with
  the number of missing cells and observed-vs-imputed mean/sd for both
  the deep-generative and classical-baseline completions.

## Details

Only the biomarker (response) columns named on the left-hand side of
`long_sub_fixed` are imputed, and they must be numeric (continuous): an
ordinal biomarker cannot be imputed with this function, and gives an
error. Every covariate referenced on the right-hand side of
`long_sub_fixed`/`long_sub_random` is assumed fully observed (an error
is raised if any covariate has missing values). The deep generative
methods model `p(biomarkers | covariates, missingness mask)` row by row
(one row per subject-visit; the covariates include the visit time
`time_variable` when `long_sub_fixed` uses it), and use
`long_sub_random` only for the id variable; `"mixed"` uses the full
random-effects structure.

Rows are aligned across biomarkers by `(id, time_variable)`: when
`data_fit_all` is a list of different data.frames with different visit
schedules per biomarker (rather than one data.frame shared by every
biomarker), the union of every `(id, time_variable)` pair seen in any of
them is used, so a visit at which one biomarker was measured but another
was not still gets a value imputed for the missing one.

Because the deep generative methods need enough data to fit reliably;
the `diagnostics` element of the return value also reports a classical
single regression imputation (each biomarker regressed on the
covariates, fit on its observed rows) as a lightweight baseline, so the
two can be compared on a given dataset before trusting the deep model's
completions.

## References

Mattei, P.-A. and Frellsen, J. (2019). MIWAE: Deep Generative Modelling
and Imputation of Incomplete Data Sets. *Proceedings of the 36th
International Conference on Machine Learning*, PMLR 97:4413-4423.

Ho, J., Jain, A., and Abbeel, P. (2020). Denoising Diffusion
Probabilistic Models. *Advances in Neural Information Processing
Systems*, 33:6840-6851.

Tashiro, Y., Song, J., Song, Y., and Ermon, S. (2021). CSDI: Conditional
Score-based Diffusion Models for Probabilistic Time Series Imputation.
*Advances in Neural Information Processing Systems*, 34:24804-24816.

## Examples

``` r
data(pbc3)
data_fit_all <- pbc3[pbc3$status3 == 1, ]

# pbc3 has no missing serBilir/albumin values; delete some at random
# (MCAR, a special case of MAR) to have something to impute
set.seed(1)
n <- nrow(data_fit_all)
data_fit_all$serBilir[sample.int(n, floor(0.1 * n))] <- NA
data_fit_all$albumin[sample.int(n, floor(0.1 * n))] <- NA

long_sub_fixed <- list(
  "long1" = serBilir ~ year + age + sex + years,
  "long2" = albumin ~ year + age + sex + years)
long_sub_random <- list("long1" = ~ year | id, "long2" = ~ year | id)

# multiple imputation from the mixed model, pooled with Rubin's rules
imputed <- imputeLongitudinal(data_fit_all, long_sub_fixed, long_sub_random,
                              time_variable = "year", impute = "multiple", seed = 1)
fits <- lapply(imputed$data_fit_all_list, longitudinalSub,
               long_sub_fixed = long_sub_fixed, long_sub_random = long_sub_random)
pooled <- poolLongitudinalSub(fits)

# \donttest{
# a deep generative method instead (needs torch)
if (requireNamespace("torch", quietly = TRUE) && torch::torch_is_installed()) {
  imputed_miwae <- imputeLongitudinal(data_fit_all, long_sub_fixed, long_sub_random,
                                      time_variable = "year", method = "miwae",
                                      epochs = 50, seed = 1)
}
# }
```
