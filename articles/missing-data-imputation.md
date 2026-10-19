# Missing data imputation

Follow-up is rarely complete: a visit happens but one lab is not drawn,
or a biomarker is only measured for some patients.
[`imputeLongitudinal()`](https://wenhaoli18.github.io/BJM/reference/imputeLongitudinal.md)
fills in the missing biomarker values before the longitudinal sub-models
are fit, and
[`poolLongitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/poolLongitudinalSub.md)
combines the fits to several completions. This article starts with the
question that comes first: **do you need to impute at all?** Often you
do not.

| Step | Function |
|----|----|
| Fill in missing biomarker values (one or several completions) | [`imputeLongitudinal()`](https://wenhaoli18.github.io/BJM/reference/imputeLongitudinal.md) |
| Fit each completion | [`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md) |
| Combine the fits (Rubin’s rules) into one, also for prediction | [`poolLongitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/poolLongitudinalSub.md) |

``` r

library(BJM)
data(pbc3)
```

## What BJM does without imputation

Without imputation,
[`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
fits each biomarker’s linear mixed model by maximum likelihood on every
visit at which *that* biomarker was measured. Under the usual assumption
that values are **missing at random** (whether a value is missing
depends only on what is observed), this is already a valid analysis: the
mixed model uses each patient’s observed values, and nothing is gained
by filling in the others with values drawn from the same model.

Only one part of the fit uses complete cases: the joint covariance of
the random effects across biomarkers, `Sigma_fit`, which every
prediction uses, is estimated from the patients with at least one value
of *every* biomarker. A patient who never had some biomarker measured is
left out of it. This is what imputation can recover.

## Do you need to impute? A simulation

To find out, we took the `pbc3` patients with an event, deleted
biomarker values at random, and compared each method’s predictions with
those of the fit to the complete data. Each method’s fits to the
completed data sets were pooled with
[`poolLongitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/poolLongitudinalSub.md)
(multiple imputation, 5 completions). The prediction error is the mean
absolute difference, over 245 patients event-free at year 3, between the
predicted 2-year risks and those of the complete-data fit; averages over
20 random deletions. In the complete data, the albumin random-intercept
variance is 0.072 and the albumin `years` (event-time) coefficient
0.078.

| Missing | Method | Prediction error | Albumin intercept variance | Albumin `years` coefficient |
|----|----|----|----|----|
| 15% of visits miss each biomarker | no imputation | **0.0056** | 0.071 | 0.078 |
|  | `"mixed"` | **0.0056** | 0.070 | 0.078 |
|  | `"miwae"`, `history = TRUE` | 0.0062 | 0.067 | 0.078 |
|  | `"diffusion"`, `history = TRUE` | 0.0077 | 0.062 | 0.076 |
|  | `"miwae"` | 0.0165 | 0.051 | 0.077 |
| 30% of patients never had albumin | no imputation | **0.0108** | 0.072 | 0.078 |
|  | `"mixed"` | 0.0114 | 0.073 | 0.078 |
|  | `"miwae"` | 0.0116 | 0.053 | 0.078 |
|  | `"miwae"`, `history = TRUE` | 0.0141 | 0.062 | 0.074 |
|  | `"diffusion"`, `history = TRUE` | 0.0157 | 0.060 | 0.067 |
| 50% of patients never had albumin | no imputation | 0.0167 | 0.070 | 0.077 |
|  | `"mixed"` | **0.0159** | 0.074 | 0.079 |
|  | `"miwae"` | 0.0192 | 0.042 | 0.078 |
|  | `"diffusion"`, `history = TRUE` | 0.0227 | 0.047 | 0.061 |
|  | `"miwae"`, `history = TRUE` | 0.0235 | 0.052 | 0.073 |

What this shows:

- **With visits missing values at random, do not impute.** No method
  beat the analysis without imputation; the mixed-model imputation tied
  it.
- **With many patients missing a biomarker entirely**, multiple
  imputation with `method = "mixed"` (the default) helped a little at
  50% (closer to the complete data in 15 of the 20 deletions), and was
  about even at 30% (closer in 8 of 20).
- **`"mixed"` is the only method that left the model unbiased.** The
  deep generative methods shrank the variation between patients (the
  albumin intercept variance), and with `history = TRUE` they also
  biased the event-time coefficient, the link between the biomarker and
  the event that every BJM prediction rests on.

Imputation is also worth it when values are missing for reasons recorded
in variables that are *not* in the model (auxiliary variables): add them
to the formulas used for imputation, then fit the model without them.

## Imputing with the mixed model

To illustrate, we delete albumin for 30% of the patients:

``` r

fixed <- list(serBilir = serBilir ~ year + age + sex + years + years:year,
              albumin  = albumin  ~ year + age + sex + years + years:year)
random <- list(serBilir = ~ year | id, albumin = ~ year | id)
complete <- pbc3[pbc3$status3 == 1, ]

set.seed(1)
never <- sample(unique(complete$id), round(0.3 * length(unique(complete$id))))
whole <- complete
whole$albumin[whole$id %in% never] <- NA
```

[`imputeLongitudinal()`](https://wenhaoli18.github.io/BJM/reference/imputeLongitudinal.md)
takes the same data and formulas as
[`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md),
plus the visit-time variable. With the default `method = "mixed"`, it
fits that same multivariate mixed model and draws each patient’s missing
values from their distribution given everything observed for the
patient: all visits, and the other biomarkers through `Sigma_fit`. For
each completion the model’s parameters are first drawn from their
approximate posterior, so the completions also carry the uncertainty
about the fit. It needs no extra package and takes a second or two.

``` r

imputed <- imputeLongitudinal(whole, fixed, random, time_variable = "year",
                              impute = "multiple", n_imputations = 5, seed = 1)
names(imputed)
#> [1] "data_fit_all"      "data_fit_all_list" "diagnostics"
```

`data_fit_all_list` holds the five completed data sets, each a list with
one data frame per biomarker, ready for
[`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
with the same formulas. Fit each, and pool the fits:

``` r

fits <- lapply(imputed$data_fit_all_list, function(d) longitudinalSub(d, fixed, random))
pooled <- poolLongitudinalSub(fits)
round(pooled$pooled$albumin, 4)
#>             estimate std_error      df statistic p_value    riv    fmi
#> (Intercept)   3.6317    0.2169 17.8601   16.7407  0.0000 0.8583 0.5135
#> year         -0.2065    0.0285 57.6136   -7.2350  0.0000 0.3309 0.2734
#> age          -0.0065    0.0027 58.9846   -2.4365  0.0179 0.2397 0.2194
#> sex          -0.2135    0.0976 13.2395   -2.1885  0.0471 1.0112 0.5640
#> years         0.0745    0.0111 33.5488    6.7271  0.0000 0.4103 0.3297
#> year:years    0.0084    0.0037 66.9882    2.2413  0.0283 0.2965 0.2507
```

`fmi`, the fraction of missing information, is the share of a
coefficient’s variance that is due to the missing values.
`pooled$long_fit_all` is a fit for prediction: the pooled fixed effects,
and the averages of the fits’ residual variances and `Sigma_fit`.
Compared with the complete-data fit and the fit without imputation:

``` r

fit_complete <- longitudinalSub(complete, fixed, random)
fit_none <- longitudinalSub(whole, fixed, random)
summarise <- function(f) c(albumin_intercept_var = f$Sigma_fit[3, 3],
                           albumin_years = nlme::fixef(f$lfit[[2]])[["years"]],
                           cor_intercepts = cov2cor(f$Sigma_fit)[1, 3])
round(rbind(complete_data = summarise(fit_complete), no_imputation = summarise(fit_none),
            mixed_imputation = summarise(pooled$long_fit_all)), 3)
#>                  albumin_intercept_var albumin_years cor_intercepts
#> complete_data                    0.072         0.078         -0.227
#> no_imputation                    0.083         0.078         -0.258
#> mixed_imputation                 0.088         0.074         -0.315
```

A single deletion like this one is noisy: here neither fit is clearly
closer to the complete data than the other. The simulation above, which
averages over 20 deletions, is what the advice rests on.

`diagnostics` compares, for each biomarker, the observed values with the
imputed ones (the means of the draws), and with a classical baseline
(each biomarker regressed on the covariates):

``` r

d <- imputed$diagnostics$albumin
round(rbind(mixed = unlist(d$deep_generative), baseline = unlist(d$classical_baseline)), 3)
#>          n_observed n_imputed observed_mean observed_sd imputed_mean imputed_sd
#> mixed           628       244         3.238       0.534        3.179      0.372
#> baseline        628       244         3.238       0.534        3.207      0.290
```

### Predicting with an imputed fit

The pooled fit is used for prediction like any other. Only the data the
models are fit on are imputed; a patient’s own history is used as
observed, missing values and all
([`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
drops a missing value and keeps the rest of the visit):

``` r

survival_fit <- survivalSub(pbc3[!duplicated(pbc3$id), ], Surv(years, status3) ~ age + sex, NULL)
patient <- pbc3[pbc3$id == 2 & pbc3$year <= 5, ]
patient$years <- NA
patient$albumin[2] <- NA
predictRisk(list(patient, patient), pooled$long_fit_all, survival_fit,
            prediction_time = 5, horizon = 2, time_variable = "year",
            survival_variable_all = list(), survival_trans_function = list())
#> 
#> =================================================================
#>  Dynamic Prediction - Event Risk
#> -----------------------------------------------------------------
#>   Competing risks   : No
#>   Subjects          : 1
#> -----------------------------------------------------------------
#> 
#>  Subject Risk Prob
#>        2    0.1531
#> 
#> =================================================================
```

### One completion or several

`impute = "single"` (the default) returns one completed data set,
`data_fit_all`, each missing value the mean of the `n_imputations`
draws. A model fit to it treats the imputed values as observed: its
standard errors are too small, and its variances, including `Sigma_fit`,
are too small the more is imputed. Prefer `impute = "multiple"` with
[`poolLongitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/poolLongitudinalSub.md)
unless very little is missing.

## Deep generative alternatives

`method = "miwae"` (a variational autoencoder for incomplete data;
Mattei and Frellsen, 2019) and `method = "diffusion"` (a conditional
diffusion model as in CSDI; Tashiro et al., 2021) learn the joint
distribution of the biomarkers and covariates with a neural network
instead. They complete each visit from that visit’s covariates and
observed biomarkers, which is why they shrink the variation between
patients; `history = TRUE` adds a summary of the patient’s other visits
to the inputs. They are flexible, but as the simulation shows, on BJM’s
longitudinal data they did not do better than the mixed model, and could
bias the model. They need the **torch** package and its backend
(`install.packages("torch")`, then
[`torch::install_torch()`](https://torch.mlverse.org/docs/reference/install_torch.html)).

``` r

torch_available <- requireNamespace("torch", quietly = TRUE) &&
  isTRUE(tryCatch(torch::torch_is_installed(), error = function(e) FALSE))
```

``` r

imputed_miwae <- imputeLongitudinal(whole, fixed, random, time_variable = "year",
                                    method = "miwae", impute = "multiple", seed = 1)
fits_miwae <- lapply(imputed_miwae$data_fit_all_list, function(d) longitudinalSub(d, fixed, random))
round(summarise(poolLongitudinalSub(fits_miwae)$long_fit_all), 3)
#> albumin_intercept_var         albumin_years        cor_intercepts 
#>                 0.063                 0.081                -0.270
```

## Things to watch for

- **Missing at random.** Every method assumes that whether a value is
  missing depends only on what is observed. A lab skipped because the
  patient was too ill breaks this, and no imputation fixes it.
- **Only biomarkers are imputed.** The covariates of the formulas and
  the visit times must be complete, and the biomarkers continuous.
- **Check on your own data.** Delete some observed values at random,
  impute them, and compare with what was deleted and with the fit
  without imputation, as in the simulation above.
- **Reproducibility.** `seed` fixes the random numbers of the call (and
  torch’s, for the deep methods, whose results can still differ slightly
  between computers).

## References

Mattei, P.-A. and Frellsen, J. (2019). MIWAE: Deep generative modelling
and imputation of incomplete data sets. *Proceedings of the 36th
International Conference on Machine Learning*, PMLR 97:4413–4423.

Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys*.
Wiley.

Tashiro, Y., Song, J., Song, Y. and Ermon, S. (2021). CSDI: Conditional
score-based diffusion models for probabilistic time series imputation.
*Advances in Neural Information Processing Systems*, 34:24804–24816.
