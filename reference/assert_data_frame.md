# Assert that an object is a non-empty data.frame

Shared input-validation helper: raises a clear error instead of letting
a malformed argument fail deep inside model-fitting code with a cryptic
message.

## Usage

``` r
assert_data_frame(x, arg_name)
```
