# Train the MIWAE network with a full-batch importance-weighted ELBO

Internal helper for
[`imputeLongitudinal`](https://wenhaoli18.github.io/BJM/reference/imputeLongitudinal.md).

## Usage

``` r
train_miwae(
  net,
  X_cov,
  X_resp,
  X_mask,
  epochs,
  K,
  learning_rate = 0.001,
  train_inputs = NULL
)
```

## Arguments

- train_inputs:

  Optional function from
  [`history_dropout_inputs()`](https://wenhaoli18.github.io/BJM/reference/history_dropout_inputs.md)
  giving the inputs with fresh random swaps of the history features each
  epoch; `NULL` uses `X_cov` as is.
