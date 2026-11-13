# Does the survival time enter the event-type model only as a main effect?

Used by
[`conditionalDT()`](https://wenhaoli18.github.io/BJM/reference/conditionalDT.md):
when the survival time appears in `form_conditional_cr` only as a plain
main effect (`status ~ time + ...`), or not at all, its effect at each
point of the integration grid is its coefficient times the grid value.
Otherwise (a transformation such as `pmin(time, 5)` or `log(time)`, or
an interaction) the model matrix has to be rebuilt at every grid point.

## Usage

``` r
time_enters_linearly(tt, survival_variable)
```

## Arguments

- tt:

  The terms of the event-type model (right-hand side).

- survival_variable:

  Name of the survival-time variable.

## Value

`TRUE` or `FALSE`.
