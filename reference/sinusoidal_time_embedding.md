# Sinusoidal diffusion-timestep embedding

Internal helper for
[`imputeLongitudinal`](https://wenhaoli18.github.io/BJM/reference/imputeLongitudinal.md)'s
`method = "diffusion"` backend: the standard Transformer/DDPM sinusoidal
position embedding, applied to a (possibly per-row-varying) integer
timestep instead of a sequence position. The frequency basis is computed
once in plain R (not with
[`torch::torch_arange()`](https://torch.mlverse.org/docs/reference/torch_arange.html),
whose inclusive/exclusive endpoint convention is a frequent source of
off-by-one bugs) and then multiplied against the timestep tensor, so the
output width is always exactly `dim` regardless of parity.

## Usage

``` r
sinusoidal_time_embedding(t, dim)
```

## Arguments

- t:

  A 1-D torch float tensor of timesteps, one per row.

- dim:

  Output embedding width.
