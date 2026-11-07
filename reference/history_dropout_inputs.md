# Random swaps of the history features while training

With `history = TRUE`, a subject never measured for a biomarker has
history features (count 0, a stand-in level) that the network almost
never sees while training, since every training row with that biomarker
observed comes from a subject who was measured. Faced with them at
imputation it extrapolated: in simulations it took the stand-in for the
subject's own level, imputed with too little spread, and biased the
fitted coefficients (the event-time coefficient by about 40%). Training
therefore swaps, at random, a row's history features for one biomarker
with their `"unseen"` version (see
[`history_features()`](https://wenhaoli18.github.io/BJM/reference/history_features.md)),
keeping the row's own values as targets, so the network learns to impute
that biomarker from the covariates and the other biomarkers alone. Each
biomarker's rate is the share of rows from subjects never measured for
it, at least 10%, and 0 (no swaps, as before) when every subject has
some value of it.

## Usage

``` r
history_dropout_inputs(X_cov, X_unseen, columns, rate)
```

## Arguments

- X_cov:

  Network inputs (standardized), a torch tensor.

- X_unseen:

  The same inputs with every biomarker's history replaced by its
  `"unseen"` version, standardized identically.

- columns:

  List with, per biomarker, the columns of its history features in
  `X_cov`.

- rate:

  Per-biomarker swap probability.

## Value

A function of a row index tensor (or `NULL` for every row) returning the
inputs with fresh random swaps, or `NULL` when every rate is 0.
