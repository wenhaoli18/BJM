# Check that a longitudinal sub-model uses the event time

Check that a longitudinal sub-model uses the event time

## Usage

``` r
assert_uses_event_time(formula_i, i, survival_variable, survival_variable_all)
```

## Arguments

- formula_i:

  The biomarker's fixed-effects formula.

## Value

Invisibly `TRUE`; an error naming biomarker `i` otherwise.
