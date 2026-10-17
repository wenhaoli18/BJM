# Build the conditional diffusion denoising network

Internal helper for
[`imputeLongitudinal`](https://wenhaoli18.github.io/BJM/reference/imputeLongitudinal.md)'s
`method = "diffusion"` backend. A single MLP predicts the noise in the
target entries of a row from: the noisy target entries (zero elsewhere),
the conditioning entries' values (zero elsewhere), the conditioning
mask, the covariates and an embedding of the diffusion step – a tabular
version of CSDI's conditional denoiser (Tashiro et al., 2021).

## Usage

``` r
build_diffusion_module(n_cov, M, hidden_units, time_emb_dim)
```

## Arguments

- n_cov:

  Number of (standardized) covariate columns.

- M:

  Number of biomarkers.

- hidden_units:

  Hidden-layer sizes of the MLP body.

- time_emb_dim:

  Width of the timestep embedding.

## Value

An instantiated
[`torch::nn_module`](https://torch.mlverse.org/docs/reference/nn_module.html).
