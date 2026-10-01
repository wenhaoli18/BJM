# Impute missing longitudinal biomarker values with a deep generative model

[`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
fits each biomarker's mixed model on the complete cases for that
biomarker, and then keeps only the subjects who have at least one
non-missing observation of *every* biomarker (see
[`longitudinalSub`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)).
With interrupted/irregular follow-up this throws away information: a
subject missing just one of several biomarkers at a visit – or missing a
biomarker's measurements entirely – is dropped from every biomarker's
fit, not just the one it is missing.

`imputeLongitudinal()` is an optional preprocessing step that fills
these gaps *before*
[`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
runs, assuming the missingness is at random (MAR) given the observed
covariates and biomarkers. Two deep generative backends are available,
selected with `method`:

- `"miwae"` (default) fits a MIWAE (Missing data Importance-Weighted
  AutoEncoder; Mattei & Frellsen, 2019) – a variational autoencoder
  whose training objective only scores reconstruction on the observed
  entries – jointly across all supplied biomarkers, and uses
  self-normalized importance resampling from the fitted decoder to draw
  one or more plausible completions of the missing cells.

- `"diffusion"` fits a conditional denoising diffusion probabilistic
  model (DDPM; Ho, Jain & Abbeel, 2020) over the same per-row biomarker
  vector, trained with a masked denoising score-matching loss that –
  like the MIWAE objective above – is only evaluated on observed entries
  (the same masking principle used by Ouyang et al.'s MissDiff, 2023,
  for training diffusion models on incomplete tabular data). Missing
  cells are then filled by a RePaint-style (Lugmayr et al., 2022)
  reverse-diffusion sampler: at every denoising step, observed
  dimensions are re-noised directly from their true value and only the
  unobserved dimensions are drawn from the learned reverse process, so
  the final sample matches the data exactly wherever it was actually
  observed.

Both backends draw one or more plausible completions of the missing
cells, and the completed data.frame(s) can be passed directly to
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
  method = c("miwae", "diffusion"),
  impute = c("single", "multiple"),
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

  Number of full-batch training epochs.

- importance_samples:

  Number of importance samples (`K`) used both in the training objective
  and at imputation time. Only used when `method = "miwae"`.

- diffusion_steps:

  Number of forward/reverse diffusion timesteps. Only used when
  `method = "diffusion"`; more steps generally improve sample quality at
  proportionally higher training/sampling cost.

- method:

  Which deep generative backend to fit: `"miwae"` (default) or
  `"diffusion"`. See Details.

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
  of imputed cells. Note also that each row (subject-visit) is imputed
  from that row's covariates and observed biomarkers only, not from the
  same subject's other visits, so imputed values do not carry a
  subject's own level or trend; this too pulls the estimated
  between-subject (random-effects) variation towards zero.

- seed:

  Optional integer seed for reproducibility. It seeds both torch and R's
  random number generator for the duration of the call; R's random
  number state from before the call is restored afterwards.

## Value

A named list with elements:

- data_fit_all:

  A completed version of `data_fit_all`, in the same list/data.frame
  shape, with every missing biomarker cell filled. Ready to pass
  straight to
  [`longitudinalSub`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md).

- data_fit_all_list:

  When `impute = "multiple"`, a list of `n_imputations` independently
  completed versions of `data_fit_all`; `NULL` when `impute = "single"`.

- diagnostics:

  A named list, one element per biomarker, each with the number of
  missing cells and observed-vs-imputed mean/sd for both the
  deep-generative and classical-baseline completions.

## Details

Only the biomarker (response) columns named on the left-hand side of
`long_sub_fixed` are imputed; every covariate referenced on the
right-hand side of `long_sub_fixed`/`long_sub_random` is assumed fully
observed (an error is raised if any covariate has missing values). Each
biomarker's own within-subject serial correlation continues to be
modeled downstream by
[`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) inside
[`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md),
exactly as today; this function only models
`p(biomarkers | covariates, missingness mask)` row by row (one row per
subject-visit), so it does not need to know the random-effects structure
in `long_sub_random` beyond the id variable it names.

Rows are aligned across biomarkers by `(id, time_variable)`: when
`data_fit_all` is a list of different data.frames with different visit
schedules per biomarker (rather than one data.frame shared by every
biomarker), the union of every `(id, time_variable)` pair seen in any of
them is used, so a visit at which one biomarker was measured but another
was not still gets a value imputed for the missing one.

Because this is a deep generative model, it needs enough data to fit
reliably; the `diagnostics` element of the return value also reports a
classical single regression imputation (each biomarker regressed on the
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

Lugmayr, A., Danelljan, M., Romero, A., Yu, F., Timofte, R., and Van
Gool, L. (2022). RePaint: Inpainting using Denoising Diffusion
Probabilistic Models. *Proceedings of the IEEE/CVF Conference on
Computer Vision and Pattern Recognition*, 11461-11471.

Ouyang, Y., Xie, L., Li, C., and Cheng, G. (2023). MissDiff: Training
Diffusion Models on Tabular Data with Missing Values. *ICML 2023
Workshop on Structured Probabilistic Inference & Generative Modeling*.

## Examples

``` r
# \donttest{
if (requireNamespace("torch", quietly = TRUE)) {
  data(pbc3)
  data_fit_all <- pbc3[pbc3$status3 == 1, ]

  # pbc3's complete-case subset has no missing serBilir/albumin values;
  # artificially delete some at random (MCAR, a special case of MAR) to
  # simulate the interrupted-follow-up gaps imputeLongitudinal() targets
  set.seed(1)
  n <- nrow(data_fit_all)
  data_fit_all$serBilir[sample.int(n, floor(0.1 * n))] <- NA
  data_fit_all$albumin[sample.int(n, floor(0.1 * n))] <- NA

  long_sub_fixed <- list(
    "long1" = serBilir ~ year + age + sex + years,
    "long2" = albumin ~ year + age + sex + years)
  long_sub_random <- list("long1" = ~ year | id, "long2" = ~ year | id)

  imputed <- imputeLongitudinal(data_fit_all, long_sub_fixed,
                                 long_sub_random, time_variable = "year",
                                 epochs = 50, seed = 1)
  long_fit_all <- longitudinalSub(imputed$data_fit_all, long_sub_fixed,
                                   long_sub_random)

  # the diffusion backend fills the same gaps via a conditional DDPM
  # instead of a VAE -- same call shape, just a different `method`
  imputed_diffusion <- imputeLongitudinal(data_fit_all, long_sub_fixed,
                                           long_sub_random, time_variable = "year",
                                           method = "diffusion", epochs = 50,
                                           diffusion_steps = 30, seed = 1)
}
#> Error: Lantern is not loaded. Please use `install_torch()` to install additional dependencies.
# }
```
