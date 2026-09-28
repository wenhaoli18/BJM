# Starting values for "auto" bandcount doubling

The built-in starting point
[`auto_tune_bandcount()`](https://liwh0904.github.io/BJM/reference/auto_tune_bandcount.md)
doubles from for each bandcount argument. `bandcount1`/
`bandcount2`/`bandcount3` used to default to fixed numbers (`10`, `40`,
and `300` respectively, across
[`predictRisk()`](https://liwh0904.github.io/BJM/reference/predictRisk.md)/[`dynamicPredictionBio()`](https://liwh0904.github.io/BJM/reference/dynamicPredictionBio.md));
those same numbers are reused here as the starting point for
auto-tuning, so that the first call
[`auto_tune_bandcount()`](https://liwh0904.github.io/BJM/reference/auto_tune_bandcount.md)
makes matches what a caller relying on the old fixed defaults would have
gotten. This cannot instead be read off `formals(predict_fun)`, because
that default is now the literal string `"auto"` itself.

## Usage

``` r
bandcount_auto_start
```

## Format

An object of class `numeric` of length 3.
