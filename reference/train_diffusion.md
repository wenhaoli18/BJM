# Train the diffusion network with a masked denoising score-matching loss

Internal helper for
[`imputeLongitudinal`](https://wenhaoli18.github.io/BJM/reference/imputeLongitudinal.md)'s
`method = "diffusion"` backend. At every epoch, each row is noised to an
independently-sampled random timestep, the network predicts the noise
that was added, and the squared-error loss is averaged only over
*observed* entries (`X_mask == 1`) – missing entries have no real `x0`
to noise in the first place, so scoring them would just train the
network to fit whatever placeholder (zero) value they were filled with.
This mirrors
[`train_miwae`](https://wenhaoli18.github.io/BJM/reference/train_miwae.md)'s
observed-entries-only reconstruction term, and matches the masked
training objective used by MissDiff (Ouyang et al., 2023) for diffusion
models on incomplete tabular data.

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
  learning_rate = 0.001
)
```
