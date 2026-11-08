# Survival sub-model that borrows a long-term tail from an external cohort

**Experimental.** Predicting beyond the last follow-up of the study
needs the baseline hazard there, which
[`survivalSub()`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md)
can only extrapolate (linearly in the cumulative hazard). When an
external cohort with longer follow-up is available – an earlier
generation of the same device, a registry, an older trial – its baseline
hazard can supply the shape of the tail. `survivalSubExternal()` fits
one Cox model to both cohorts with a generation indicator (1 for the
study, 0 for the external cohort): \$\$h(t \mid x, g) = h_0(t)
\exp(\beta^\top x + \gamma g),\$\$ so the study's hazard beyond its own
follow-up is the external one scaled by \\e^\gamma\\. Only the survival
sub-model uses the external cohort; the longitudinal sub-model is fit on
the study as usual, and the external cohort needs no biomarkers.

The key assumption is proportional hazards between the cohorts, also
beyond the study's follow-up, where it cannot be checked. Within the
overlap it is tested with
[`cox.zph`](https://rdrr.io/pkg/survival/man/cox.zph.html) (stored in
`$external$ph_test`, with a warning when \\p \< 0.05\\); the test has
little power with short study follow-up, so a non-significant result is
weak reassurance, and a sensitivity analysis on the tail is advisable.

Predictions are for the study population: when the prediction data have
no `generation` column it is taken as 1.

## Usage

``` r
survivalSubExternal(
  data_internal,
  data_external,
  form_marginal_surv,
  form_conditional_cr = NULL,
  generation = "generation"
)
```

## Arguments

- data_internal, data_external:

  One row per subject, for the study and for the external cohort. Both
  need the variables of `form_marginal_surv` (and of
  `form_conditional_cr`, if given).

- form_marginal_surv, form_conditional_cr:

  As in
  [`survivalSub`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md),
  without the generation term, which is added.

- generation:

  Name of the generation indicator column to create (it must not already
  exist).

## Value

A `survivalSub.BJM` fit, usable wherever
[`survivalSub()`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md)'s
is (e.g.
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md),
`fitImputedBJM(survival_fit_all = )`), with an `external` element: the
generation variable, the cohort sizes and follow-up, the hazard ratio of
the study against the external cohort and the PH test.

## Examples

``` r
# \donttest{
data(pbc3)
surv <- pbc3[!duplicated(pbc3$id), ]
# pretend the study had 4 years of follow-up and the rest of pbc3 is external
id <- as.numeric(as.character(surv$id))
study <- surv[id %% 2 == 0, ]
study$status3[study$years > 4] <- 0
study$years <- pmin(study$years, 4)
external <- surv[id %% 2 == 1, ]
fit <- survivalSubExternal(study, external, Surv(years, status3) ~ age + sex, NULL)
fit$external
#> $generation
#> [1] "generation"
#> 
#> $n_internal
#> [1] 156
#> 
#> $n_external
#> [1] 156
#> 
#> $max_followup_internal
#> [1] 4
#> 
#> $max_followup_external
#> [1] 14.30566
#> 
#> $hazard_ratio
#> [1] 0.6608038
#> 
#> $ph_test
#>     chisq        df         p 
#> 0.5637548 1.0000000 0.4527513 
#> 
# }
```
