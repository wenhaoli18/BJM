# Nonparametric estimate of the event-time distribution from interval data

Self-consistency (EM) estimate in the spirit of Turnbull (1976) for
interval-censored event times \\(L, R\]\\, left-truncated at \\V\\ (each
subject is only in the sample because \\T \> V\\), with an optional
known event type per detected event (the competing-risks estimate of
Hudgens, Satten & Longini, 2001, Biometrics 57:74). Mass sits on the
cells between consecutive distinct time points (all \\L\\, \\R\\, \\V\\
and `extra`); within an observed interval it is not identified, and the
EM's choice there is what the estimate reports.

## Usage

``` r
npmle_interval(
  L,
  R,
  V,
  cause = NULL,
  extra = NULL,
  tol = 1e-09,
  max_iter = 5000
)
```

## Arguments

- L, R, V:

  Interval and truncation time per subject (`R = Inf` when no event was
  detected).

- cause:

  Event type per subject (1..K, `NA` without an event), or `NULL`.

- extra:

  Additional time points to put cell boundaries at.

- tol, max_iter:

  EM convergence settings.

## Value

A list with `edges` and `mass` (cells x event types, the last cell
running to `Inf`).
