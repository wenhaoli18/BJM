# Reshape recurrent-event data into gap-time episodes for the backward joint model

**Experimental.** Prepares recurrent-event data for the backward joint
model by a gap-time (renewal) reformulation: it splits each subject's
event history into *episodes* – the interval from one restart (the
previous event, or the discharge from the stay it began) to the next
event, with the clock reset to `0` at the start of each – and re-clocks
any accompanying longitudinal measurements onto the same per-episode
timeline. The result plugs into the existing
[`survivalSub()`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md)/
[`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)/[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
pipeline, unchanged, to predict for a subject currently at risk in an
open episode whether their *next* event (within a horizon) is a
recurrence or a terminal event;
[`predictRecurrentCount`](https://wenhaoli18.github.io/BJM/reference/predictRecurrentCount.md)
builds on the same data to predict the *number* of recurrences in a
window. Unlike the single-event and competing-risks methods the rest of
the package rests on (Shen and Li, 2021; Li et al., 2026 – see
`citation("BJM")`), this extension is unpublished, and has been checked
only in simulation studies, not on real data (see the "Recurrent events"
article on the package website).

## Usage

``` r
buildRecurrentEpisodes(
  data_events,
  data_long = NULL,
  id_variable = "id",
  time_variable = "time",
  status_variable = "status",
  time_variable_long = NULL,
  episode_variable = "episode",
  gap_variable = "gap",
  time_in_episode_variable = "time_in_episode",
  event_any_variable = "event_any",
  event_type_variable = "event_type",
  restart_variable = NULL,
  lookback = 0,
  episode_id_variable = "episode_id"
)
```

## Arguments

- data_events:

  One row per subject per event-or-censoring occurrence, in calendar
  time (i.e. not yet re-clocked to gap time – this function does that).
  Must contain `id_variable`, `time_variable` (numeric, strictly
  increasing within each subject), and `status_variable` (`0` =
  censored/end of follow-up, `1` = recurrent event, `2` = terminal event
  – the same 0/1/2 convention
  [pbc3](https://wenhaoli18.github.io/BJM/reference/pbc3.md)'s `status5`
  already uses). Any other columns (baseline or event-time-varying
  covariates) are carried through to `episode_data`, taken from the same
  row as the event they end. Subjects are matched across `data_events`
  and `data_long` by the printed value of their id, so an id that is a
  number in one table and a factor or character in the other still
  matches; `id` keeps the type it has in `data_events`.

- data_long:

  Optional repeated-measurement data to re-clock onto the same
  per-episode timeline, in the same format
  [`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)/
  [`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
  accept elsewhere in this package: a single bare `data.frame`, or a
  list of `data.frame`s (one per longitudinal outcome). Must contain
  `id_variable` and `time_variable_long`. Rows whose
  `time_variable_long` falls outside every episode – before `0` (beyond
  `lookback`), in a hospital stay `lookback` does not reach, or after
  the subject's process has already ended (per `data_events`) – are
  dropped, with a message reporting how many. A row whose
  `time_variable_long` falls exactly on the shared boundary between two
  episodes (e.g. a biomarker drawn at the same visit an event is
  recorded) is genuinely ambiguous – it is at once "the last measurement
  of the ending episode" and "the first (`time_in_episode = 0`)
  measurement of the next one" – so it is kept *twice*, once per
  episode, rather than silently assigned to only one (which would
  systematically discard every fresh-episode baseline measurement that
  happens to coincide with the previous episode's event time, and bias
  the fitted backward-covariate effect); a message says how many rows
  this applies to.

- id_variable, time_variable, status_variable:

  Column names in `data_events` (defaults `"id"`, `"time"`, `"status"`).

- time_variable_long:

  Column name of the measurement time in `data_long`; defaults to
  `time_variable` (i.e. the same name as in `data_events`) when `NULL`.

- episode_variable, gap_variable, time_in_episode_variable,
  event_any_variable, event_type_variable, episode_id_variable:

  Column names to create in the output (defaults `"episode"`, `"gap"`,
  `"time_in_episode"`, `"event_any"`, `"event_type"`, `"episode_id"`).
  `episode_id` is `"<id>:<episode>"` (a factor), one value per episode:
  the grouping variable to use for the biomarker model's random effects
  (see Details). Change these if they collide with an existing column of
  `data_events`/`data_long` – collisions are an error, not a silent
  overwrite.

- restart_variable:

  Optional name of a numeric column of `data_events`: for each
  recurrence (`status == 1`) row, the time the process restarts after it
  (for hospitalizations, the discharge ending that stay; must not be
  earlier than the row's `time_variable`). The next episode starts there
  instead of at the event time. `NULL` (default) restarts at the event
  time itself. A missing value on a recurrence row that starts another
  episode is taken as no stay, with a warning. Ignored on rows that end
  the process. Not carried into `episode_data`.

- lookback:

  A single non-negative number (`Inf` allowed): how far before an
  episode's start its measurements are still kept, as negative
  `time_in_episode`, but never from before the previous event. Default
  `0` keeps only measurements from the start on.

## Value

A list with elements:

- episode_data:

  One row per subject per episode: `id_variable`, `episode_variable`
  (episode number, starting at 1), `gap_variable` (that episode's gap
  time), `event_any_variable` (`0`/`1`: did the episode end in an event
  of either type, i.e. is `gap` an observed event time rather than an
  administrative censoring time), `event_type_variable` (`0` =
  recurrence, `1` = terminal event, `NA` when `event_any == 0`), and
  every other column of `data_events`.

- long_data:

  `NULL` if `data_long` was not supplied; otherwise `data_long`, in the
  same bare-`data.frame`-or-list shape it was supplied in, with
  `episode_variable`, `gap_variable`, `time_in_episode_variable`,
  `event_any_variable`, and `event_type_variable` columns added (and
  out-of-range rows dropped, as described above). Any other column of
  `episode_data` (baseline or event-time-varying covariates, e.g.
  `age`/`sex`) is also carried onto matching rows here, the way every
  row of [pbc3](https://wenhaoli18.github.io/BJM/reference/pbc3.md)'s
  long-format table already repeats `age`/`sex`/`years`/`status3` – but
  only into a column `data_long` does not already have of that name.

## Details

**What counts as an episode.** Episode 1 runs from time `0` to a
subject's first row in `data_events`; episode `k > 1` runs from the end
of episode `k - 1` to that subject's `k`-th row. A subject contributes
no further episodes after a row with `status %in% c(0, 2)`
(administrative censoring or a terminal event); any later rows for that
`id` are ignored, with a message saying how many. `status == 1` rows
(recurrences) always start a new episode.

**Feeding the result into
[`survivalSub()`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md).**
Fit the gap-time recurrence hazard and the recurrence-vs-terminal-event
type model exactly as
[`survivalSub()`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md)
already fits a marginal hazard and a competing-risks type model, using
`episode_data` in place of one-row-per-subject data,
`Surv(gap, event_any) ~ covariates` in place of `form_marginal_surv`,
and `event_type ~ gap + covariates` in place of `form_conditional_cr`
(`event_type` is `0` for a recurrence, `1` for a terminal event – so
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)'s
`risk_prob_1`/`risk_prob_2` come out as the recurrence/terminal risk
respectively). Administratively-censored episodes get `event_any = 0`
and `event_type = NA`, exactly like `status3 = 0`/`status4 = NA` for a
censored subject in
[`pbc3`](https://wenhaoli18.github.io/BJM/reference/pbc3.md). The
episodes of one subject are then treated as independent: the
coefficients are still consistent for the gap-time model, but the
standard errors that [`summary()`](https://rdrr.io/r/base/summary.html)
reports ignore the correlation between a subject's episodes and are too
small.

**Feeding the result into
[`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)/[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md).**
`gap` plays the same dual role `years` plays for
[pbc3](https://wenhaoli18.github.io/BJM/reference/pbc3.md): it is both
the [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) time in
`episode_data` and the backward covariate swept over candidate values
inside
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md),
so it must appear in `long_sub_fixed` (e.g.
`biomarker ~ time_in_episode + gap + gap * time_in_episode + covariates`),
and the training data passed to
[`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
must be filtered to `event_any == 1` rows of `long_data` – only
completed episodes have an observed `gap` – exactly as the package's own
examples filter `pbc3` to `status3 == 1` before fitting
[`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md).
`time_in_episode` (elapsed time since the current episode began) plays
the role `year` plays for `pbc3`. To predict from "now" for a subject
currently at risk in an open episode, add a trailing `status = 0` row to
`data_events` at the current data-cutoff time, and call
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
with `prediction_time` set to that subject's elapsed `time_in_episode`
at the cutoff.

**Random effects: group by episode, not by patient.** Fit the biomarker
model with the random effects grouped by `episode_id`
(`random = ~ time_in_episode | episode_id`), so that each episode is its
own "subject", as it already is in the survival sub-model. A patient's
biomarker level generally differs from one episode to the next (the
underlying trajectory keeps moving), so a single random effect per
patient is the wrong unit; in a simulation in which the biomarker's
current level drives the risk, grouping by patient gave a clearly worse
AUC than grouping by episode.
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
and
[`predictRecurrentCount()`](https://wenhaoli18.github.io/BJM/reference/predictRecurrentCount.md)
then take one `episode_id` per prediction, so give each prediction its
own `id`.

**The open episode's `gap` is `NA`.** For an episode that has not ended
(`event_any == 0`) – typically the one a subject is currently at risk in
– the `gap` column of `long_data` is `NA`: its final length is exactly
what is being predicted.
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
reads a missing survival-time variable as "event time not yet known,
still at risk", so that episode's rows can be passed straight to it,
with no placeholder to set or override (`episode_data` keeps the real
elapsed time of such an episode, which is what the Cox fit needs as a
censoring time). Training rows are unaffected, since
[`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
is fit on `event_any == 1` rows only. Note that
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
also expects the outcome columns it was fit with (the time and
event-type variables) to exist in the prediction data; `long_data`
already carries `event_type` (`NA` for an open episode).

**Episodes that start at discharge (`restart_variable`, `lookback`).**
By default an episode starts at the previous event (time `0` for the
first). For hospitalization data that is not the at-risk clock: after a
readmission the patient is in hospital, not at risk, until discharge.
Give `data_events` a numeric column – passed as `restart_variable` –
holding, for each recurrence row, the time the process restarts (the
discharge ending that stay); the next episode then starts there and its
gap is measured from it. Time `0` is the first episode's start (e.g. the
index discharge). Biomarkers measured during the stay, or before the
index discharge, precede the episode they predict: `lookback` (default
`0`, none) keeps measurements up to that long before an episode's start
– never reaching back past the previous event, so a stay's own labs are
not mixed into the earlier at-risk period – and gives them *negative*
`time_in_episode`. Use `Inf` to keep the whole stay (and, for the first
episode, everything before time `0`). Negative times are fine for
[`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
and
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md).

## Examples

``` r
## A minimal example: reshaping alone (no model fitting).
## Subject 1: two recurrences, then a terminal event; subject 2: censored.
data_events <- data.frame(
  id     = c(1, 1, 1, 2),
  time   = c(1.2, 2.5, 4.0, 3.1),
  status = c(1,   1,   2,   0),
  age    = c(58,  58,  58,  61)
)
episodes <- buildRecurrentEpisodes(data_events)
episodes$episode_data
#>   id episode episode_id gap event_any event_type age
#> 1  1       1        1:1 1.2         1          0  58
#> 2  1       2        1:2 1.3         1          0  58
#> 3  1       3        1:3 1.5         1          1  58
#> 4  2       1        2:1 3.1         0         NA  61

## Hospitalizations: each episode starts at DISCHARGE, not at the readmission.
## `time` = readmission (or death/censoring) day; `discharge` = day the stay
## it started ends. Labs drawn in the stay (and before the index discharge at
## day 0) are kept with negative time_in_episode via lookback = Inf.
events <- data.frame(id = 1, time = c(20, 60, 90), status = c(1, 1, 0),
                     discharge = c(26, 66, NA))
labs <- data.frame(id = 1, day = c(-3, 0, 5, 21, 26, 40, 61, 66, 80))
hf <- buildRecurrentEpisodes(events, labs, time_variable_long = "day",
                             restart_variable = "discharge", lookback = Inf)
hf$episode_data   # gaps measured from each discharge: 20, 34, 24
#>   id episode episode_id gap event_any event_type
#> 1  1       1        1:1  20         1          0
#> 2  1       2        1:2  34         1          0
#> 3  1       3        1:3  24         0         NA
hf$long_data[, c("day", "episode", "time_in_episode", "gap")]
#>   day episode time_in_episode gap
#> 1  -3       1              -3  20
#> 2   0       1               0  20
#> 3   5       1               5  20
#> 4  21       2              -5  34
#> 5  26       2               0  34
#> 6  40       2              14  34
#> 7  61       3              -5  NA
#> 8  66       3               0  NA
#> 9  80       3              14  NA

# \donttest{
## A full pipeline on simulated recurrent-event + longitudinal data:
## build episodes -> survivalSub() -> longitudinalSub() -> predictRisk() and
## predictRecurrentCount().
set.seed(1)
n <- 60
data_events <- do.call(rbind, lapply(seq_len(n), function(id) {
  k <- sample(1:3, 1)
  gaps <- rexp(k, rate = 1)
  data.frame(
    id = id, time = cumsum(gaps),
    status = c(rep(1, k - 1), sample(c(0, 2), 1, prob = c(0.6, 0.4))),
    age = round(50 + stats::rnorm(1, sd = 8)), sex = id %% 2
  )
}))

## subject 9999 is still at risk right now, 0.7 time units into their
## 2nd episode (1 prior recurrence, no resolution yet) -- this is who
## we will predict for below.
data_events <- rbind(data_events, data.frame(
  id = 9999, time = c(0.8, 1.6), status = c(1, 0), age = 55, sex = 1
))

data_long <- do.call(rbind, lapply(c(seq_len(n), 9999), function(id) {
  maxt <- max(data_events$time[data_events$id == id]) + 0.2
  tt <- seq(0, maxt, by = 0.25)
  data.frame(id = id, year = tt, biomarker = rnorm(length(tt)))
}))

built <- buildRecurrentEpisodes(data_events, data_long, time_variable_long = "year")
#> buildRecurrentEpisodes(): dropped 50 row(s) of `data_long` whose year fell outside every modeled episode (before time 0 or the `lookback` window, inside a hospital stay not covered by `lookback`, or after the subject's recurrent-event process had already ended).

## gap-time recurrence hazard + recurrence-vs-terminal-event type model
survival_fit_recur <- survivalSub(
  built$episode_data,
  Surv(gap, event_any) ~ age + sex + episode,
  event_type ~ gap + age + sex + episode
)

## longitudinal sub-model, trained on completed episodes only, with the
## random effects grouped by episode (see Details)
long_fit_recur <- longitudinalSub(
  built$long_data[built$long_data$event_any == 1, ],
  biomarker ~ time_in_episode + age + sex + gap + gap * time_in_episode,
  ~ time_in_episode | episode_id
)

## subject 9999's open episode: its gap is NA, which predictRisk() reads as
## "still at risk" -- nothing to override
open <- built$long_data[built$long_data$id == 9999 & built$long_data$event_any == 0, ]
landmark <- max(open$time_in_episode)

## the risk that their *next* event, within 1 time unit, is a recurrence /
## a terminal event
risk <- predictRisk(
  open, long_fit_recur, survival_fit_recur,
  prediction_time = landmark, horizon = 1, time_variable = "time_in_episode",
  survival_variable_all = list(), survival_trans_function = list(),
  bandcount1 = 10, bandcount2 = 10
)
risk$risk_prob_1 # P(next event within horizon is a recurrence)
#>    9999:2 
#> 0.2978591 
risk$risk_prob_2 # P(next event within horizon is the terminal event)
#>    9999:2 
#> 0.2763965 

## the number of recurrences over the next 3 time units
predictRecurrentCount(
  open, long_fit_recur, survival_fit_recur,
  prediction_time = landmark, horizon = 3, time_variable = "time_in_episode",
  n_sim = 100, seed = 1
)
#> Predicted recurrences in (0.7, 3.7], 1 patient(s), 100 simulated paths each
#> mean expected count 0.470; P(>=1) 0.380; P(>=2) 0.080; P(death) 0.890
#> 
#>  episode_id mean_count p_count_0 p_count_1 p_count_2 p_count_3plus p_ge1 p_ge2
#>      9999:2       0.47      0.62       0.3      0.07          0.01  0.38  0.08
#>  p_death
#>     0.89
# }
```
