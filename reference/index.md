# Package index

## Fit sub-models

Fit the survival and longitudinal sub-models that dynamic prediction is
built on.

- [`survivalSub()`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md)
  : Fit the survival sub-model
- [`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
  : Fit a multivariate longitudinal sub-model

## Dynamic prediction

Predict event risk and future biomarker values from the fitted
sub-models.

- [`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
  : Dynamic prediction function for future event risk
- [`predictLongitudinal()`](https://wenhaoli18.github.io/BJM/reference/predictLongitudinal.md)
  : Dynamic prediction function for future longitudinal outcomes
- [`checkBandcountConvergence()`](https://wenhaoli18.github.io/BJM/reference/checkBandcountConvergence.md)
  : Check whether bandcount1/bandcount2/bandcount3 are large enough
- [`survivalTrans()`](https://wenhaoli18.github.io/BJM/reference/survivalTrans.md)
  : Build a survival-time transform basis from cut points

## Missing data imputation

Impute missing longitudinal biomarker values and pool sub-model fits
across imputations.

- [`imputeLongitudinal()`](https://wenhaoli18.github.io/BJM/reference/imputeLongitudinal.md)
  : Impute missing longitudinal biomarker values with a deep generative
  model
- [`poolLongitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/poolLongitudinalSub.md)
  : Pool longitudinal sub-model fits across multiple imputations with
  Rubin's rules

## Plots

Visualize predictions, fitted sub-models, observed longitudinal
trajectories, and observed event outcomes.

- [`predictPlot()`](https://wenhaoli18.github.io/BJM/reference/predictPlot.md)
  : Plot risk and future biomarker predictions across a horizon sweep
- [`riskPlot()`](https://wenhaoli18.github.io/BJM/reference/riskPlot.md)
  : Plot predicted risk across a sweep of landmark times
- [`cmtPlot()`](https://wenhaoli18.github.io/BJM/reference/cmtPlot.md) :
  Plot conditional mean trajectories (CMT)
- [`spaghettiPlot()`](https://wenhaoli18.github.io/BJM/reference/spaghettiPlot.md)
  : Plot individual longitudinal trajectories (spaghetti plot)
- [`cifPlot()`](https://wenhaoli18.github.io/BJM/reference/cifPlot.md) :
  Plot cumulative incidence functions
- [`performancePlot()`](https://wenhaoli18.github.io/BJM/reference/performancePlot.md)
  : Plot predictive performance across landmark times
- [`calibrationPlot()`](https://wenhaoli18.github.io/BJM/reference/calibrationPlot.md)
  : Plot calibration of dynamic risk predictions
- [`plot(`*`<longitudinalSub.BJM>`*`)`](https://wenhaoli18.github.io/BJM/reference/plot.longitudinalSub.BJM.md)
  : Plot a fitted longitudinal sub-model
- [`plot(`*`<survivalSub.BJM>`*`)`](https://wenhaoli18.github.io/BJM/reference/plot.survivalSub.BJM.md)
  : Plot a fitted survival sub-model
- [`plot(`*`<predictLongitudinal.BJM>`*`)`](https://wenhaoli18.github.io/BJM/reference/plot.predictLongitudinal.BJM.md)
  : Plot predicted biomarker distributions

## Data

- [`pbc3`](https://wenhaoli18.github.io/BJM/reference/pbc3.md) : Mayo
  Clinic primary biliary cirrhosis data used as example code
