# Recode ordinal biomarker responses onto the training categories

Shared helper for `predictRisk`, `dynamicPredictionBio`, and
`dynamicPredictionBioAll`. The copula densities turn an ordinal response
into its category code with
[`as.numeric()`](https://rdrr.io/r/base/numeric.html) and index the
fitted thresholds with it, and the candidate categories of an ordinal
`bio_i` are read off [`levels()`](https://rdrr.io/r/base/levels.html).
Both used to come from the prediction data itself, so any difference
from the training categories – unused levels dropped (e.g. by
[`droplevels()`](https://rdrr.io/r/base/droplevels.html) or subsetting
and re-creating the factor), a different level order, or a character
column – silently bracketed the latent score between the wrong
thresholds (in one check, a risk of 0.53 instead of 0.10). Each ordinal
response is now rebuilt as an ordered factor on the categories the
`clmm()` fit was estimated with, matching by label; a value that is not
one of those categories is an error. Missing values stay missing
([`drop_missing_longitudinal()`](https://wenhaoli18.github.io/BJM/reference/drop_missing_longitudinal.md)
removes them afterwards).

## Usage

``` r
align_ordinal_levels(data_predict_all, long_fit_all)
```

## Value

`data_predict_all`, with each ordinal response recoded.
