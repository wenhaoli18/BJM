# Package index

## Fit sub-models

Fit the survival and longitudinal sub-models that dynamic prediction is
built on.

- [`survivalSub()`](https://liwh0904.github.io/BJM/reference/survivalSub.md)
  : Fit the survival sub-model
- [`longitudinalSub()`](https://liwh0904.github.io/BJM/reference/longitudinalSub.md)
  : Fit a multivariate longitudinal sub-model

## Dynamic prediction

Predict event risk and future biomarker values from the fitted
sub-models.

- [`predictRisk()`](https://liwh0904.github.io/BJM/reference/predictRisk.md)
  : Dynamic prediction function for future event risk
- [`predictLongitudinal()`](https://liwh0904.github.io/BJM/reference/predictLongitudinal.md)
  : Dynamic prediction function for future longitudinal outcomes
- [`checkBandcountConvergence()`](https://liwh0904.github.io/BJM/reference/checkBandcountConvergence.md)
  : Check whether bandcount1/bandcount2/bandcount3 are large enough
- [`survivalTrans()`](https://liwh0904.github.io/BJM/reference/survivalTrans.md)
  : Build a survival-time transform basis from cut points

## Missing data imputation

Impute missing longitudinal biomarker values and pool sub-model fits
across imputations.

- [`imputeLongitudinal()`](https://liwh0904.github.io/BJM/reference/imputeLongitudinal.md)
  : Impute missing longitudinal biomarker values with a deep generative
  model
- [`poolLongitudinalSub()`](https://liwh0904.github.io/BJM/reference/poolLongitudinalSub.md)
  : Pool longitudinal sub-model fits across multiple imputations with
  Rubin's rules

## Plots

Visualize predictions and observed longitudinal trajectories.

- [`predictPlot()`](https://liwh0904.github.io/BJM/reference/predictPlot.md)
  : Plot risk and future biomarker predictions across a horizon sweep
- [`riskPlot()`](https://liwh0904.github.io/BJM/reference/riskPlot.md) :
  Plot predicted risk across a sweep of landmark times
- [`cmtPlot()`](https://liwh0904.github.io/BJM/reference/cmtPlot.md) :
  Plot conditional mean trajectories (CMT)

## Data

- [`pbc3`](https://liwh0904.github.io/BJM/reference/pbc3.md) : Mayo
  Clinic primary biliary cirrhosis data used as example code
