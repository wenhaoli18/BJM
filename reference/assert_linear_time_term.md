# Check that the survival time enters the event-type model linearly

Used by
[`survivalSub()`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md):
at prediction time the event-type model `form_conditional_cr` is
evaluated at every point of the survival-time integration grid by
multiplying the survival time's coefficient by the grid value. That is
only correct when the survival time appears as a plain main effect
(`status ~ time + ...`); a transformation (`log(time)`, `I(time^2)`) or
an interaction (`time:age`) would silently be evaluated at the wrong
value, so it is rejected here. Leaving the survival time out entirely is
allowed (event type then does not depend on event time).

## Usage

``` r
assert_linear_time_term(form_conditional_cr, survival_variable)
```

## Arguments

- form_conditional_cr:

  The event-type formula.

- survival_variable:

  Name of the survival-time variable.
