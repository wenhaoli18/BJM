# Baseline cumulative hazard spliced at the study's last follow-up

For `survivalSubExternal(splice = TRUE)`: the Breslow estimate from the
study alone, at the joint fit's coefficients (with the study's
generation term), up to `tau`; beyond it, the joint fit's baseline
increments. Both are on the scale of the joint fit's uncentered linear
predictor, so
[`survival_lp()`](https://wenhaoli18.github.io/BJM/reference/survival_accessors.md)
applies unchanged.

## Usage

``` r
splice_baseline(fit, data_internal, y_int, tau)
```
