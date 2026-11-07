# Predict the number of recurrent events in a window, by rolling episodes forward

**Experimental.** For a patient at the start of (or partway through) an
episode of a gap-time recurrent-event model fit with
[`buildRecurrentEpisodes`](https://wenhaoli18.github.io/BJM/reference/buildRecurrentEpisodes.md),
predicts the distribution of the *number* of recurrences in
`(prediction_time, prediction_time + horizon]`, together with the
probability of dying in it, by Monte Carlo: every patient is simulated
forward `n_sim` times, episode after episode, until the window closes or
the patient dies.
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
answers "is the *next* event a recurrence or a terminal event, and
when"; this answers "how many recurrences".

## Usage

``` r
predictRecurrentCount(
  data_predict,
  long_fit_all,
  survival_fit_all,
  prediction_time,
  horizon,
  time_variable,
  survival_variable_all = list(),
  survival_trans_function = list(),
  episode_variable = "episode",
  next_covariates = NULL,
  stay_lengths = 0,
  n_sim = 200,
  max_stages = 10,
  bandcount2 = 100,
  seed = NULL
)
```

## Arguments

- data_predict:

  The rows of the patients' *open* episode, as
  [`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
  takes them (for example the `long_data` rows of
  [`buildRecurrentEpisodes()`](https://wenhaoli18.github.io/BJM/reference/buildRecurrentEpisodes.md)
  for each patient's last episode, whose `gap` is `NA`): one
  `data.frame` used for every biomarker, or a list of data frames, one
  per biomarker, as
  [`buildRecurrentEpisodes()`](https://wenhaoli18.github.io/BJM/reference/buildRecurrentEpisodes.md)
  returns for a list `data_long`; one patient per `id`: more than one
  episode for the same patient is an error. A patient whose `gap` is not
  `NA` and lies below `prediction_time` counts as no longer at risk, and
  is an error too rather than a silent zero. The patient's last row (in
  the first data frame that has the patient) supplies the covariates of
  later episodes.

- long_fit_all, survival_fit_all:

  Fits from
  [`longitudinalSub`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
  and
  [`survivalSub`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md),
  on episode data.

- prediction_time:

  The time (on the episode clock, as in
  [`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md))
  the prediction is made at: `0` at discharge.

- horizon:

  Length of the window after `prediction_time`, in the same time unit.

- time_variable:

  Name of the time variable of the longitudinal sub-model (e.g.
  `"time_in_episode"`).

- survival_variable_all, survival_trans_function:

  As in
  [`predictRisk`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md).

- episode_variable:

  Name of the episode-counter covariate (default `"episode"`); used only
  by the default `next_covariates`, and ignored if the data has no such
  column.

- next_covariates:

  Optional `function(row, k)` returning the one-row covariate data frame
  of the episode `k` episodes after the open one (`k = 1, 2, ...`),
  given the patient's last row of `data_predict` as `row`. The default
  adds `k` to `episode_variable`, but not beyond the largest episode
  number in the data the survival sub-model was fit on: past it the
  fitted effect of the counter is an extrapolation (a hazard that grows
  with each earlier recurrence would be multiplied a hundredfold or more
  by episode 30), and the simulated paths would run on through ever
  shorter episodes.

- stay_lengths:

  A numeric vector of non-negative stay lengths (the time between a
  recurrence and the discharge that restarts the process), resampled
  with replacement. `0` (default): none.

- n_sim:

  Paths simulated per patient.

- max_stages:

  Most episodes drawn per path (including the first); paths still
  running at the cap are counted as they stand, with a warning.

- bandcount2:

  Event-time grid size passed to
  [`simulateTrajectory`](https://wenhaoli18.github.io/BJM/reference/simulateTrajectory.md).

- seed:

  Optional integer seed; the caller's random number stream is left
  unchanged.

## Value

An object of class `"predictRecurrentCount.BJM"`, a list with

- summary:

  One row per patient: `mean_count` (expected number of recurrences in
  the window), `p_count_0`, `p_count_1`, `p_count_2`, `p_count_3plus`,
  `p_ge1`, `p_ge2`, and `p_death` (probability of dying in the window,
  outside hospital: a death during a stay is not part of the model, see
  [`buildRecurrentEpisodes`](https://wenhaoli18.github.io/BJM/reference/buildRecurrentEpisodes.md)).

- counts:

  The raw draws: a patients x `n_sim` matrix of recurrence counts.

- died:

  A patients x `n_sim` logical matrix: died in the window.

- horizon, prediction_time, n_sim:

  As supplied.

- stages:

  Most episodes drawn on any path.

- capped:

  `TRUE` if some path was still running when `max_stages` stopped it (a
  warning, giving their share, is also given): its later recurrences are
  not counted, so the counts are slightly low.

- capped_fraction:

  The share of all simulated paths stopped by `max_stages` (`0` if
  none).

## Details

Each simulated path goes as follows.

1.  The *current* (open) episode's event time and type are drawn with
    [`simulateTrajectory`](https://wenhaoli18.github.io/BJM/reference/simulateTrajectory.md),
    from their posterior given the patient's biomarker history in that
    episode – the same distribution
    [`predictRisk`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
    integrates, so the first step uses the biomarkers exactly as
    [`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
    does.

2.  If the event is a terminal event inside the window, the path ends
    (the patient dies). If it is a recurrence inside the window, the
    count goes up by one; a stay length is then drawn (from
    `stay_lengths`) and the next episode starts at discharge. If the
    event (or the discharge) falls after the window closes, the path
    ends.

3.  Each *later* episode is drawn the same way, but with no biomarker
    history – the patient has not been observed in it yet – so from the
    model's prior given the patient's covariates, with the episode
    counter advanced by `next_covariates`. Repeat until no path is still
    running, or `max_stages` episodes have been drawn.

**Limitation.** Only the first step uses the patient's biomarkers. Later
episodes depend on the patient only through their covariates (and the
episode counter), so any patient-level tendency to recur that is not in
those covariates – a frailty – is not carried forward. The patients who
reach later episodes are then mostly the frailer ones, the fitted effect
of the episode counter absorbs that selection, and rolled forward it
makes the expected counts too extreme: too high for the patients
predicted to recur most, too low for those predicted to recur least (in
simulations with a strong frailty the calibration slope of the expected
count was 0.7-0.8). Putting the patient's observed gap history among the
survival sub-model's covariates – an indicator of the first episode and
the mean log length of the earlier gaps, say – recovers part of it
(slope 0.83-0.9). Keep such covariates at their observed values for the
later episodes: the default `next_covariates` changes only the episode
counter. See the recurrent-events article.

The stay length is not part of the backward joint model, so it has to be
supplied. The default, `0`, means a recurrence takes no time (a tumour
recurrence, say); for hospitalizations pass observed stay lengths (they
are resampled with replacement).

## See also

[`buildRecurrentEpisodes`](https://wenhaoli18.github.io/BJM/reference/buildRecurrentEpisodes.md),
[`predictRisk`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md),
[`simulateTrajectory`](https://wenhaoli18.github.io/BJM/reference/simulateTrajectory.md)

## Examples

``` r
## Simulated recurrent events and a biomarker. Subject 9999 is still at
## risk in their 2nd episode, and is who we predict for.
set.seed(1)
n <- 60
data_events <- do.call(rbind, lapply(seq_len(n), function(id) {
  k <- sample(1:3, 1)
  data.frame(
    id = id, time = cumsum(rexp(k, rate = 1)),
    status = c(rep(1, k - 1), sample(c(0, 2), 1, prob = c(0.6, 0.4))),
    age = round(50 + rnorm(1, sd = 8)), sex = id %% 2
  )
}))
data_events <- rbind(data_events, data.frame(
  id = 9999, time = c(0.8, 1.6), status = c(1, 0), age = 55, sex = 1
))
data_long <- do.call(rbind, lapply(c(seq_len(n), 9999), function(id) {
  tt <- seq(0, max(data_events$time[data_events$id == id]) + 0.2, by = 0.25)
  data.frame(id = id, year = tt, biomarker = rnorm(length(tt)))
}))
built <- buildRecurrentEpisodes(data_events, data_long, time_variable_long = "year")
#> buildRecurrentEpisodes(): dropped 50 row(s) of `data_long` whose year fell outside every modeled episode (before time 0 or the `lookback` window, inside a hospital stay not covered by `lookback`, or after the subject's recurrent-event process had already ended).

survival_fit_recur <- survivalSub(
  built$episode_data,
  Surv(gap, event_any) ~ age + sex + episode,
  event_type ~ gap + age + sex + episode
)
long_fit_recur <- longitudinalSub(
  built$long_data[built$long_data$event_any == 1, ],
  biomarker ~ time_in_episode + age + sex + gap + gap * time_in_episode,
  ~ time_in_episode | episode_id
)

## The number of recurrences over the next 3 time units
open <- built$long_data[built$long_data$id == 9999 & built$long_data$event_any == 0, ]
counts <- predictRecurrentCount(
  open, long_fit_recur, survival_fit_recur,
  prediction_time = max(open$time_in_episode), horizon = 3,
  time_variable = "time_in_episode", n_sim = 50, bandcount2 = 30, seed = 1
)
counts
#> Predicted recurrences in (0.7, 3.7], 1 patient(s), 50 simulated paths each
#> mean expected count 0.480; P(>=1) 0.380; P(>=2) 0.080; P(death) 0.880
#> 
#>  episode_id mean_count p_count_0 p_count_1 p_count_2 p_count_3plus p_ge1 p_ge2
#>      9999:2       0.48      0.62       0.3      0.06          0.02  0.38  0.08
#>  p_death
#>     0.88
```
