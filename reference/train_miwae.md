# Train the MIWAE network with a full-batch importance-weighted ELBO

Internal helper for
[`imputeLongitudinal`](https://liwh0904.github.io/BJM/reference/imputeLongitudinal.md).

## Usage

``` r
train_miwae(net, X_cov, X_resp, X_mask, epochs, K, learning_rate = 0.001)
```
