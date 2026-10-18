# Interactive demo

This interactive demo runs BJM’s dynamic prediction for a single patient
in the browser, with no R needed. Pick a patient from the `pbc3` data,
or upload one patient’s follow-up visits as a CSV file. Then choose a
prediction time and a window and press **Predict**. The interface is
available in English and Chinese.

If the frame does not load, [open the demo in a new
tab](https://wenhaoli18.shinyapps.io/bjm-demo/). On the free hosting
plan the app goes to sleep when nobody is using it, so the first visit
can take up to a minute to start.

## What the demo shows

The app fits the two sub-models once, when it starts, on `pbc3`: a Cox
model for death (`status3`) with age and sex, and linear mixed models
for three biomarkers (log serum bilirubin, albumin and log prothrombin
time). For the selected patient it then uses only the visits up to the
prediction time $`t`$ and shows:

| Panel | BJM function |
|----|----|
| Cumulative risk of the event in $`(t, t + h]`$ for $`h`$ up to the window | [`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md) |
| Predicted mean and 95% interval of each biomarker at the end of the window | [`predictLongitudinal()`](https://wenhaoli18.github.io/BJM/reference/predictLongitudinal.md) |
| Observed history and predicted future trajectory of each biomarker | [`predictLongitudinal()`](https://wenhaoli18.github.io/BJM/reference/predictLongitudinal.md) |

The biomarker predictions are conditional on the patient still being
event-free at each future time. This is the quantity the backward joint
model predicts jointly with the event risk.

## Uploading a patient

The CSV needs one row per visit and the columns `year`, `age`, `sex`,
`serBilir`, `albumin` and `prothrombin`, on the same scales as `pbc3`.
The app has a template you can download. The event time is not needed:
BJM integrates over it.

## Running the app locally

The app’s source is in the
[`shiny/bjm-demo`](https://github.com/wenhaoli18/BJM/tree/master/shiny/bjm-demo)
folder of the repository. To run it on your own computer, clone the
repository and run:

``` r

install.packages(c("shiny", "bslib"))
shiny::runApp("shiny/bjm-demo")
```

The demo is for illustrating the method only. Its models are fit on the
PBC data set and must not be used for clinical decisions.
