# Build a linear variance-preserving diffusion (DDPM) noise schedule

Internal helper for
[`imputeLongitudinal`](https://liwh0904.github.io/BJM/reference/imputeLongitudinal.md)'s
`method = "diffusion"` backend. Returns the standard Ho et al. (2020)
linear `beta` schedule and the derived `alpha`/ `alpha_bar` sequences,
as plain numeric vectors (not torch tensors) so
[`diffusion_impute`](https://liwh0904.github.io/BJM/reference/diffusion_impute.md)
can index them by a plain R integer timestep inside its reverse-sampling
loop without repeated tensor\<-\>R round trips.

## Usage

``` r
diffusion_beta_schedule(diffusion_steps, beta_start = 1e-04, beta_end = 0.02)
```
