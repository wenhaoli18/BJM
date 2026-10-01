# Mode of a density tabulated on a grid

Helper for
[`compute_bio_marker_step()`](https://wenhaoli18.github.io/BJM/reference/compute_bio_marker_step.md):
the grid point where `density` peaks, refined (when `refine = TRUE`) by
fitting a parabola through the log-density at that point and its two
neighbours. For a smooth, near-normal predictive density this recovers
the mode to well within one grid step, instead of snapping to the grid:
the unrefined mode moved in whole grid steps as `bandcount3` changed,
which for a predicted value near 0 is a large relative change, so
`"auto"` tuning of `bandcount3` reported non-convergence even on a fine
grid. No refinement at the grid's ends, for an ordinal biomarker (whose
grid is its categories), or if a neighbour's density is 0.

## Usage

``` r
density_mode(Y_all, density, refine = TRUE)
```

## Arguments

- Y_all:

  The (equally spaced, for continuous biomarkers) grid.

- density:

  The density at each grid point.

- refine:

  Whether to refine between grid points.

## Value

A single number, or `NA` if `density` has no finite maximum.
