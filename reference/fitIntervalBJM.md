# Fit the backward joint model to interval-censored event times

**Experimental.** With interval censoring nobody's event time \\T\\ is
known exactly – only \\L \< T \le R\\ – so the longitudinal sub-model
\\f(Y \mid T)\\, which has \\T\\ as a covariate, cannot be fit on the
subjects with an observed event the way it is for right-censored data.
`fitIntervalBJM()` fills in \\T\\ by stochastic EM, then by multiple
imputation:

1.  fit the marginal survival model to the intervals
    ([`survivalSub()`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md)
    with `Surv(L, R, type = "interval2")`); it does not depend on \\T\\,
    so it is fit once;

2.  start every interval-censored \\T\\ at its interval midpoint and fit
    [`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md);

3.  draw each subject's \\T\\ from \\f(T \mid L \< T \le R, Y, X)
    \propto f(Y \mid T, X) f(T \mid X)\\ on \\(L, R\]\\
    ([`imputeEventTime()`](https://wenhaoli18.github.io/BJM/reference/imputeEventTime.md)),
    refit
    [`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md),
    and repeat for `n_burnin` iterations;

4.  keep iterating, saving every `n_thin`-th draw and its
    [`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
    fit until `n_imputations` are saved.

Step 3 uses the backward model's own factorization: \\f(Y \mid T)\\ on a
grid of \\T\\ is exactly what dynamic prediction already computes
([`conditionalYT()`](https://wenhaoli18.github.io/BJM/reference/conditionalYT.md)),
so each draw costs about one
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
call per subject.

For a subject whose event was detected, only the longitudinal
measurements up to the last visit without the event (\\L\\) are used:
the one at the detecting visit \\R\\ (and any later one) is taken after
the event, which \\f(Y \mid T)\\ does not describe, and they are left
out with a message. As in the rest of the interval-censoring support,
the biomarker visits are taken to be the inspection visits, so no
measurement lies between \\L\\ and \\R\\.

By default only subjects whose event was observed (finite \\R\\) enter
the longitudinal fit, mirroring the complete-case fit used for
right-censored data, which is valid when censoring is independent of
\\(T, Y)\\ given the covariates. `include_right_censored = TRUE` also
imputes \\T \in (L, \infty)\\ for right-censored subjects and uses them;
their draws then rely on the survival model's extrapolated tail.

The saved draws come from one chain, so successive draws are correlated;
saving only every `n_thin`-th iteration makes them closer to
independent. They are drawn with the model parameters at their point
estimates, not from their posterior, so the between-imputation variance
in `pooled` can still be somewhat understated.

For prediction, `long_fit_all` combines the imputations: its fixed
effects are the Rubin's-rules pooled estimates and its residual and
random-effects variances the averages over the imputations (see
[`poolLongitudinalSub`](https://wenhaoli18.github.io/BJM/reference/poolLongitudinalSub.md)),
so a prediction does not depend on which single imputation happened to
come first. With competing risks the event-type model of
`survival_fit_all` likewise has its coefficients averaged over the
imputations.

## Usage

``` r
fitIntervalBJM(
  data_survival,
  data_fit_all,
  form_marginal_surv,
  event_time,
  long_sub_fixed,
  long_sub_random,
  time_variable,
  survival_variable_all = NULL,
  survival_trans_function = NULL,
  form_conditional_cr = NULL,
  n_burnin = 10,
  n_imputations = 5,
  n_thin = 5,
  n_grid = 40,
  baseline = c("spline", "piecewise"),
  df = NULL,
  include_right_censored = FALSE,
  seed = NULL
)
```

## Arguments

- data_survival:

  One row per subject: the interval endpoints, the survival covariates
  and the subject id.

- data_fit_all:

  As in
  [`longitudinalSub`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md):
  one data frame per biomarker (or one shared data frame), long format,
  for every subject. The `event_time` column and the columns named in
  `survival_variable_all` are (re)filled from the imputed \\T\\.

- form_marginal_surv:

  `Surv(L, R, type = "interval2") ~ covariates`.

- event_time:

  Name of the event-time column used in `long_sub_fixed`.

- long_sub_fixed, long_sub_random:

  As in
  [`longitudinalSub`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md).

- time_variable:

  Name of the visit-time column.

- survival_variable_all, survival_trans_function:

  Transformed event-time columns, as in
  [`predictRisk`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
  (`NULL` if none).

- form_conditional_cr:

  Optional competing-risks (event-type) model, as in
  [`survivalSub`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md),
  e.g. `type ~ Tev + x`. The event type (0/1, in `data_survival` and in
  every `data_fit_all` data frame, alongside the covariates of this
  formula) is taken as known once the event is detected; each subject's
  \\T\\ is then drawn from \\f(Y \mid T, D) P(D \mid T) f(T)\\, and this
  model, which uses \\T\\, is refit on every draw.

- n_burnin:

  Stochastic-EM iterations before draws are saved.

- n_imputations:

  Number of saved draws (and longitudinal fits).

- n_thin:

  Iterations between saved draws (1 saves every iteration after the
  burn-in). Larger values make the saved draws less correlated at the
  cost of `n_burnin + n_imputations * n_thin` iterations in all.

- n_grid:

  Grid points per subject interval.

- baseline, df:

  Baseline hazard of the survival model (see
  [`survivalSub`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md)).

- include_right_censored:

  See Description.

- seed:

  Optional integer seed (the global RNG state is restored).

## Value

An object of class `"fitIntervalBJM"`, a list with `survival_fit_all`
(with competing risks, its event-type model has the coefficients
averaged over the imputations, while the rest of that `glm` object, e.g.
its standard errors, is the first imputation's fit;
`survival_fit_all_list` has one fit per imputation); `long_fit_all_list`
(one
[`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
fit per imputation); `long_fit_all` (the imputations combined, see
Description, for use with
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
and friends); `pooled`
([`poolLongitudinalSub`](https://wenhaoli18.github.io/BJM/reference/poolLongitudinalSub.md)
of the fits, when `n_imputations > 1`); `imputed_T` (a subjects x
imputations matrix); and `trace` (each iteration's fixed effects, to
check the chain has settled).

## Examples

``` r
# \donttest{
# simulated study: visits every 1-3 years, an event known only to lie
# between the last visit without it (L) and the first visit with it (R)
set.seed(1)
n <- 60
x <- rnorm(n)
event_time <- (-log(runif(n)) / (0.05 * exp(0.5 * x)))^(1 / 1.5)
surv <- data.frame(id = seq_len(n), x = x, L = NA_real_, R = NA_real_)
long <- list()
for (i in seq_len(n)) {
  visits <- c(0, cumsum(runif(30, 1, 3)))
  visits <- visits[visits < 12]
  if (event_time[i] > max(visits)) {
    surv$L[i] <- max(visits)
    seen <- visits
  } else {
    k <- findInterval(event_time[i], visits)
    surv$L[i] <- visits[k]
    surv$R[i] <- visits[k + 1]
    seen <- visits[visits <= event_time[i]]
  }
  long[[i]] <- data.frame(id = i, year = seen, x = x[i],
                          y = 1 + 0.3 * seen - 0.2 * event_time[i] + 0.4 * x[i] +
                            rnorm(1, 0, 0.5) + rnorm(length(seen), 0, 0.3))
}
long <- do.call(rbind, long)

fit <- fitIntervalBJM(surv, long, Surv(L, R, type = "interval2") ~ x,
                      event_time = "Tev", long_sub_fixed = y ~ year + Tev + x,
                      long_sub_random = ~ year | id, time_variable = "year",
                      n_burnin = 2, n_imputations = 2, n_thin = 1, n_grid = 15,
                      seed = 1)
fit$pooled
#> 
#> Call:
#> poolLongitudinalSub(m = 2 completed-data fits, M = 1 longitudinal outcome)
#> 
#> =================================================================
#>  Rubin's-Rules-Pooled Longitudinal Sub-model(s)
#> -----------------------------------------------------------------
#>  [1] Outcome 1
#> 
#>                Value       SE    DF      t  p-value    
#> (Intercept)  0.85680  0.19666 36.31  4.357 0.000104 ***
#> year         0.29168  0.01648 92.54 17.696  < 2e-16 ***
#> Tev         -0.16460  0.03693 20.72 -4.458 0.000224 ***
#> x            0.38435  0.09815 45.43  3.916 0.000300 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>   Fraction of missing information (FMI), by coefficient:
#>   (Intercept) = 0.169,  year = 0.033,  Tev = 0.222,  x = 0.051 
#> 
#> =================================================================

# the fit plugs into the prediction functions; Tev is unknown for a new patient
history <- long[long$id == 1 & long$year <= 2, ]
history$Tev <- NA_real_
predictRisk(history, fit$long_fit_all, fit$survival_fit_all,
            prediction_time = 2, horizon = 2, time_variable = "year",
            survival_variable_all = NULL, survival_trans_function = NULL,
            bandcount1 = 10, bandcount2 = 20)
#> 
#> =================================================================
#>  Dynamic Prediction - Event Risk
#> -----------------------------------------------------------------
#>   Competing risks   : No
#>   Subjects          : 1
#> -----------------------------------------------------------------
#> 
#>  Subject Risk Prob Already, undetected
#>        1    0.3007               0.592
#> 
#>   Interval-censored: risks condition on being event-free at each subject's
#>  last visit; 'Already, undetected' is the probability the event happened
#>  between that visit and the prediction time.
#> 
#> =================================================================
#> 
# }
```
