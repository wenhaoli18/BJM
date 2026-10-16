# Average the event-type model over the imputations

With competing risks
[`fitIntervalBJM()`](https://wenhaoli18.github.io/BJM/reference/fitIntervalBJM.md)
refits the event-type `glm` on every imputation, since it uses the event
time. For prediction the coefficients are averaged over the imputations
(the Rubin's-rules point estimate); the rest of the returned `glm`
object is the first imputation's fit.

## Usage

``` r
pool_event_type_fit(survival_fit_all_list)
```

## Arguments

- survival_fit_all_list:

  One `survivalSub.BJM` fit per imputation.

## Value

The first fit, with the averaged event-type coefficients.
