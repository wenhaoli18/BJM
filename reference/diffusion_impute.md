# Draw conditional reverse-diffusion completions

Internal helper for
[`imputeLongitudinal`](https://wenhaoli18.github.io/BJM/reference/imputeLongitudinal.md)'s
`method = "diffusion"` backend. The missing entries of every row are the
targets and all its observed entries the conditions: starting from pure
noise on the targets, the DDPM reverse process (Ho et al., 2020) is run
with the network conditioned on the observed values at every step.
Observed entries are returned unchanged.

## Usage

``` r
diffusion_impute(
  net,
  X_cov,
  X_resp,
  X_mask,
  resp_mean,
  resp_sd,
  schedule,
  diffusion_steps,
  time_emb_dim,
  n_imputations
)
```

## Arguments

- net:

  Trained network.

- X_cov, X_resp, X_mask:

  As in
  [`train_diffusion()`](https://wenhaoli18.github.io/BJM/reference/train_diffusion.md).

- resp_mean, resp_sd:

  Means and standard deviations that undo the responses'
  standardization.

- schedule, diffusion_steps, time_emb_dim:

  As used in training.

- n_imputations:

  Number of completions to draw.

## Value

A list of `n_imputations` completed response matrices on the original
scale.
