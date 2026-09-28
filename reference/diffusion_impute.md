# Draw RePaint-style conditional reverse-diffusion completions

Internal helper for
[`imputeLongitudinal`](https://liwh0904.github.io/BJM/reference/imputeLongitudinal.md)'s
`method = "diffusion"` backend: implements the RePaint (Lugmayr et al.,
2022) inpainting sampler, adapted from image patches to a per-row
biomarker vector. Starting from pure noise, at every reverse step `t`,
observed dimensions are overwritten by re-noising the true
(standardized) observed value directly to noise level `t - 1` (the
"known" branch, sampled from the forward process posterior, not the
network), while unobserved dimensions are drawn from the network's
learned reverse step (the "unknown" branch); the two are recombined with
the mask before the next step. This guarantees the final sample matches
the data exactly at every previously-observed cell, and only genuinely
extrapolates the missing ones. Already-observed cells are left untouched
by the caller, matching
[`miwae_impute`](https://liwh0904.github.io/BJM/reference/miwae_impute.md)'s
contract.

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
