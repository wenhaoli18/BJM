# Train the conditional diffusion network

Internal helper for
[`imputeLongitudinal`](https://wenhaoli18.github.io/BJM/reference/imputeLongitudinal.md)'s
`method = "diffusion"` backend, following CSDI's self-supervised
training (Tashiro et al., 2021): in every mini-batch, each row's
observed entries are split at random into targets (each with a
probability drawn uniformly per row) and conditions; the targets are
noised to a random diffusion step and the network learns to predict that
noise from the conditions' true values. Only observed entries ever serve
as targets, so the loss never touches a missing cell. An earlier version
noised whole rows without giving the network the other observed values,
so imputations made little use of the row's other biomarkers.

## Usage

``` r
train_diffusion(
  net,
  X_cov,
  X_resp,
  X_mask,
  epochs,
  diffusion_steps,
  time_emb_dim,
  batch_size = 128,
  learning_rate = 0.001
)
```

## Arguments

- net:

  Network from
  [`build_diffusion_module()`](https://wenhaoli18.github.io/BJM/reference/build_diffusion_module.md).

- X_cov, X_resp, X_mask:

  Covariates, standardized responses (0 where missing) and the
  observation mask, as `torch` tensors.

- epochs:

  Passes over the data.

- diffusion_steps:

  Number of diffusion steps.

- time_emb_dim:

  Width of the timestep embedding.

- batch_size:

  Rows per mini-batch.

- learning_rate:

  Adam learning rate.

## Value

A list with the trained `net` and its noise `schedule`.
