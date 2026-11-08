# Many censored subjects

The backward joint model factorizes the joint distribution as
$`f(T)\,f(Y \mid T)`$. The usual fit estimates the longitudinal
sub-model $`f(Y \mid T)`$ on the subjects with an observed event, since
only they have a known $`T`$. Prediction, however, needs $`f(Y \mid T)`$
for every $`T`$ after the prediction time. When many subjects are
censored, the events are mostly early ones, so for a patient who is
doing well the event-time effect is extrapolated well beyond the event
times it was estimated from. With a term linear in $`T`$, this pulls the
predictive distribution of $`T`$ towards early events: predicted risks
are too high, and more so at later landmarks, where the patients still
at risk are the long survivors.

[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
warns (once per session) when the longitudinal fit has no more subjects
than the survival fit has events while at least 30% of the subjects are
censored.

## Using the censored subjects

A censored subject’s history says that $`T > C`$ and how the biomarkers
behaved until $`C`$, which is exactly the information about long event
times that the events-only fit leaves out.
[`fitImputedBJM()`](https://wenhaoli18.github.io/BJM/reference/fitImputedBJM.md)
uses it: it draws each censored subject’s $`T`$ from
$`f(T \mid T > C, Y, X) \propto f(Y \mid T,
X) f(T \mid X)`$ (with competing risks, $`T`$ and the event type
jointly), refits
[`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
on everyone, iterates, and pools several imputations, as
[`fitIntervalBJM()`](https://wenhaoli18.github.io/BJM/reference/fitIntervalBJM.md)
does for interval censoring.

Once the censored subjects are in the fit, the event-time effect is
estimated over long event times too, and a term linear in $`T`$ can
drift over them. A bounded effect such as $`\exp(-T/\kappa)`$ levels off
instead. Write it in place of `years` in every term:

``` r

library(BJM)
data(pbc3)
pbc3_surv <- pbc3[!duplicated(pbc3$id), ]
fixed <- function(b, kappa) {
  G <- sprintf("I(exp(-years/%g))", kappa)
  as.formula(sprintf(paste0("%s ~ year + age + sex + %2$s + year:%2$s + status4 + status4:year + ",
                            "status4:%2$s + status4:year:%2$s"), b, G))
}
biomarkers <- c("serBilir", "albumin", "prothrombin")
fit_kappa <- function(kappa) {
  fx <- lapply(biomarkers, fixed, kappa = kappa)
  rd <- rep(list(~ year | id), 3)
  names(fx) <- names(rd) <- biomarkers
  fitImputedBJM(pbc3_surv, pbc3, Surv(years, status3) ~ age + sex, fx, rd, "year",
                form_conditional_cr = status4 ~ years + age + sex, seed = 1)
}
fits <- lapply(c(2, 5, 10), fit_kappa)
```

## Choosing the form of the event-time effect

[`observedLongLogLik()`](https://wenhaoli18.github.io/BJM/reference/observedLongLogLik.md)
gives the observed-data log-likelihood of the longitudinal sub-model:
$`\log f(Y \mid T, D)`$ for a subject with an event,
$`\log \int_C^\infty \sum_d f(Y \mid t, d) P(d \mid t)\,dF(t) - \log S(C)`$
for a censored one. Fits sharing the survival sub-model can be compared
on it; higher is better.

``` r

sapply(fits, observedLongLogLik, data_survival = pbc3_surv, data_fit_all = pbc3,
       time_variable = "year")
```

On `pbc3`, $`\kappa = 5`$ has the highest log-likelihood (about
$`-1670`$, against $`-1690`$ to $`-1720`$ for 2 and 10, and $`-2010`$
for the usual linear fit on the events).

## What it changes on pbc3

Predicted over observed 3-year cumulative incidence of death (transplant
competing; Aalen-Johansen), 5-fold cross-validation, 3 biomarkers:

| Landmark (years) | 1 | 2 | 3 |
|----|----|----|----|
| Usual fit (events only, linear in `years`) | 1.08 | 1.14 | 1.32 |
| [`fitImputedBJM()`](https://wenhaoli18.github.io/BJM/reference/fitImputedBJM.md), $`\exp(-T/5)`$ | 0.99 | 1.02 | 1.07 |

The AUC is unchanged (within 0.01). Most of the correction comes from
using the censored subjects; the bounded form keeps it from depending on
how far a linear term is extrapolated.

## Beyond the last follow-up

The draws of $`T`$ beyond the last follow-up, and any prediction
reaching there, rest on the survival sub-model’s baseline hazard, which
[`survivalSub()`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md)
can only extrapolate. When an external cohort with longer follow-up is
available (an earlier generation of a device, a registry),
[`survivalSubExternal()`](https://wenhaoli18.github.io/BJM/reference/survivalSubExternal.md)
fits one Cox model to both with a generation indicator, so the study’s
hazard there is the external one scaled by a hazard ratio. It needs
proportional hazards between the cohorts, also where only the external
cohort is observed; the `cox.zph` test it runs on the overlap has little
power with short follow-up, so treat the borrowed tail as an assumption
to vary in a sensitivity analysis.

``` r

# sketch: study_surv / study_long are the study (one row per subject / long format),
# external_surv the external cohort (one row per subject, no biomarkers needed);
# fx, rd as above
surv_ext <- survivalSubExternal(study_surv, external_surv, Surv(years, status3) ~ age + sex)
surv_ext$external$ph_test
fit <- fitImputedBJM(study_surv, study_long, Surv(years, status3) ~ age + sex, fx, rd, "year",
                     survival_fit_all = surv_ext)
```

## The whole predictive distribution of $`T`$

[`predictEventTime()`](https://wenhaoli18.github.io/BJM/reference/predictEventTime.md)
returns, for each patient, the CDF of $`T`$ at several horizons and its
quantiles in one pass. The 10% quantile, for example, is a time the
patient outlives with probability 0.9 – a lower prediction bound for
durability.

``` r

history <- pbc3[pbc3$id %in% c(2, 7) & pbc3$year <= 3, ]
history$years <- NA
history$status4 <- NA
predictEventTime(history, fits[[2]]$long_fit_all, fits[[2]]$survival_fit_all,
                 prediction_time = 3, horizons = c(1, 3, 5), time_variable = "year",
                 probs = c(0.1, 0.5, 0.9))
```
