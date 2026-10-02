# Expand cumulative incidence confidence limits into a step-shaped ribbon

`geom_ribbon()` joins points with straight lines; to follow the step
curves, each jump time gets two rows, one with the previous limits and
one with the new limits.

## Usage

``` r
cif_step_ribbon(curves)
```

## Arguments

- curves:

  The curve `data.frame` built in
  [`cifPlot()`](https://wenhaoli18.github.io/BJM/reference/cifPlot.md).

## Value

A `data.frame` with columns `time`, `lower`, `upper`, `hue`, `curve`.
