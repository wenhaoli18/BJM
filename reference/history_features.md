# What a subject's other visits say about each biomarker

Inputs for the deep imputation models (`history = TRUE` in
[`imputeLongitudinal`](https://wenhaoli18.github.io/BJM/reference/imputeLongitudinal.md)):
without them each visit is completed from that visit alone, so an
imputed value does not carry the subject's own level, which shrinks the
estimated variation between subjects. For each visit and biomarker: the
mean of the subject's *other* observed values of the biomarker, how many
there are (`log1p`), and the value at the nearest other observed visit
and its distance in time (`log1p`). The visit's own value is left out,
so the model cannot simply copy it while training. A subject with no
other observed value gets the biomarker's overall mean, a count of 0 and
the largest distance.

## Usage

``` r
history_features(aligned, id, time_variable, resp_mat, mask)
```

## Value

A matrix with four columns per biomarker, one row per visit.
