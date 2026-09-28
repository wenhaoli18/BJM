# Assert that an index is a valid, in-range biomarker position

Shared input-validation helper for `bio_i`: catches an out-of-range or
non-integer value before it becomes a "subscript out of bounds" error
from indexing into `data_predict_all`/`lfit`.

## Usage

``` r
assert_index(x, max_value, arg_name, context)
```
