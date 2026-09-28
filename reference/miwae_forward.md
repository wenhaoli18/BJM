# Draw K importance-weighted latent samples and decode them

Internal helper shared by
[`train_miwae`](https://liwh0904.github.io/BJM/reference/train_miwae.md)
and
[`miwae_impute`](https://liwh0904.github.io/BJM/reference/miwae_impute.md):
encodes each row once, draws K samples from q(z\|x_o) (reusing the same
encoder mean/logvar for all K, as in the original MIWAE), decodes all of
them, and returns everything needed to compute the importance-weighted
log-likelihood terms.

## Usage

``` r
miwae_forward(net, X_cov, X_resp, X_mask, K)
```
