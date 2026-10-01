# Number of subjects the random-effects covariance is estimated from

The EM update for `Sigma_fit` averages each subject's \\E\[b b^T\]\\
over the subjects that are actually in the fit – those with at least one
complete observation of every biomarker (`yi` has one element per such
subject). It used to divide by the number of subjects in the first
biomarker's raw data instead, which underestimated `Sigma_fit` (badly,
since EM compounds it over iterations) whenever some subject was
excluded, e.g. one never measured for some biomarker.

## Usage

``` r
n_subjects(yi)
```

## Arguments

- yi:

  The per-subject list of stacked responses.

## Value

The number of subjects.
