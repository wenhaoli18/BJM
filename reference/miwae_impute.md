# Draw self-normalized-importance-resampled completions from a fitted MIWAE

Internal helper for
[`imputeLongitudinal`](https://wenhaoli18.github.io/BJM/reference/imputeLongitudinal.md):
implements the MIWAE multiple-imputation procedure (Mattei & Frellsen,
2019, Section 3.2) – draw K decoder samples per row, importance-weight
them by their (masked-to-observed-entries) likelihood, then resample one
index per completion and use that draw's decoder mean plus Gaussian
noise as the imputed value. Already-observed cells are left untouched by
the caller; this only needs to return sensible values for missing cells.

## Usage

``` r
miwae_impute(net, X_cov, X_resp, X_mask, resp_mean, resp_sd, K, n_imputations)
```
