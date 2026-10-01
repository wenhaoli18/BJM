# Seed R's RNG until the calling function returns

Calls `set.seed(seed)` and registers, in the calling function's frame,
an exit handler that puts back the `.Random.seed` that existed before
(or removes it if there was none), so a `seed` argument does not change
the caller's random number stream.

## Usage

``` r
local_r_seed(seed, frame = parent.frame())
```

## Arguments

- seed:

  Integer seed.

- frame:

  The frame whose exit restores the RNG state.
