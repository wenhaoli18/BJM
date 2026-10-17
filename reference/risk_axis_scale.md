# Biomarker axis with the risk probabilities on a secondary axis

Helper for
[`riskPlot()`](https://wenhaoli18.github.io/BJM/reference/riskPlot.md)
and
[`predictPlot()`](https://wenhaoli18.github.io/BJM/reference/predictPlot.md),
which draw risks on the biomarker's axis multiplied by `scale_prob`. The
secondary (risk) axis is labelled only within \\\[0, 1\]\\: where the
biomarker goes below 0 it used to show negative "probabilities". For an
ordinal biomarker the primary axis is labelled with its categories.

## Usage

``` r
risk_axis_scale(scale_prob, category_labels = NULL)
```

## Arguments

- scale_prob:

  Factor the risks were multiplied by.

- category_labels:

  Category labels of an ordinal biomarker (plotted as their codes
  `1:K`), or `NULL`.

## Value

A ggplot2 scale.
