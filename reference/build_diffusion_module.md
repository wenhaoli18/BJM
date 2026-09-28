# Build the diffusion denoising network

Internal helper for
[`imputeLongitudinal`](https://wenhaoli18.github.io/BJM/reference/imputeLongitudinal.md)'s
`method = "diffusion"` backend. A single MLP predicts the noise `eps`
added to the (standardized) biomarker vector, conditioned on the noisy
input itself, the observed-data mask, the fully-observed covariates, and
an embedding of the current timestep – the same "condition the
generative model on fully-observed covariates and the mask" design used
by
[`build_miwae_module`](https://wenhaoli18.github.io/BJM/reference/build_miwae_module.md).

## Usage

``` r
build_diffusion_module(n_cov, M, hidden_units, time_emb_dim)
```
