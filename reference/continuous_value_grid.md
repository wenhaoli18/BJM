# Candidate-value grid for a continuous biomarker

Helper for
[`compute_bio_marker_step()`](https://wenhaoli18.github.io/BJM/reference/compute_bio_marker_step.md):
the grid of candidate values (`Y_all`) a continuous biomarker's
predictive density is tabulated on, shared by all patients predicted.
Built in two passes:

1.  A coarse grid over a deliberately wide range – 5 spans either side
    of the patients' observed values, the span being their range but at
    least the biomarker's SD in the training data – with a step of the
    biomarker's residual SD \\\sigma\\, at least 30 and at most 100
    points. A predicted future measurement includes its measurement
    error, so every predictive density has an SD of at least \\\sigma\\:
    within half a step of its mode it is still above 88\\ `rel_tol`, so
    no patient's density falls between coarse points.

2.  `bandcount3` + 1 equally spaced points over the part of the coarse
    grid where some patient's density exceeds `rel_tol` times that
    patient's maximum, widened by \\2\sigma\\ (and at least one coarse
    step) on either side.

The wide range used to be the final grid, so most of the `bandcount3`
points (about 90\\ every density was practically 0. If the coarse pass
finds no usable density, the wide range is returned with `bandcount3` +
1 points, as before.

## Usage

``` r
continuous_value_grid(
  shared,
  bio_i,
  bio_i_name,
  long_fit_all,
  bandcount3,
  density_fun,
  rel_tol = 1e-06
)
```

## Arguments

- shared:

  Output of
  [`compute_bio_shared_step()`](https://wenhaoli18.github.io/BJM/reference/compute_bio_shared_step.md).

- bio_i:

  Index of the biomarker.

- bio_i_name:

  Its response variable name.

- long_fit_all:

  Output of
  [`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md).

- bandcount3:

  Number of intervals of the final grid.

- density_fun:

  Function of a vector of candidate values returning the density matrix
  (one row per value, one column per patient).

- rel_tol:

  Relative density below which a value is outside the grid.

## Value

A numeric vector of candidate values.
