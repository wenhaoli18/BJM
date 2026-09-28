# Build the MIWAE encoder/decoder network

Internal helper for
[`imputeLongitudinal`](https://liwh0904.github.io/BJM/reference/imputeLongitudinal.md).
The encoder input is (standardized covariates, zero-filled standardized
biomarkers, missingness mask); the decoder input is (latent draw,
standardized covariates), matching the "condition the generative model
on fully-observed covariates too" design used throughout.

## Usage

``` r
build_miwae_module(n_cov, M, hidden_units, latent_dim)
```
