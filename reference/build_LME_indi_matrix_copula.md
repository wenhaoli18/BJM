# Build one biomarker's transposed design matrix at prediction time (Gaussian-copula-aware)

Copula-aware counterpart to the
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) call
inline in
[`conditionalYT()`](https://liwh0904.github.io/BJM/reference/conditionalYT.md)/[`conditionalYDT()`](https://liwh0904.github.io/BJM/reference/conditionalYDT.md).
Continuous markers use the cached `terms`/`xlevels`/`contrasts` from the
[`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) fit, exactly as
those functions do. Ordinal markers instead reuse
[`longitudinalSubCopula()`](https://liwh0904.github.io/BJM/reference/longitudinalSubCopula.md)'s
own construction verbatim (`model.matrix(long_sub_fixed[[i]], data)`
subset to `names(lfit[[i]]$beta)`), since an
[`ordinal::clmm`](https://rdrr.io/pkg/ordinal/man/clmm.html) fit's
`$terms` does not carry the no-intercept convention `lfit[[i]]$beta`
expects.

## Usage

``` r
build_LME_indi_matrix_copula(i, data_i, lfit, long_fit_all, long_sub_fixed)
```

## Value

A matrix with one row per fixed-effect coefficient and one column per
observation (row of `data_i`).
