# Mayo Clinic primary biliary cirrhosis data used as example code

The dataset originates from the Mayo Clinic trial on primary biliary
cirrhosis (PBC) of the liver, carried out from 1974 to 1984. It includes
data from 424 PBC patients who were referred to the Mayo Clinic within
this decade and met the eligibility requirements for a randomized
placebo-controlled trial of D-penicillamine. However, only the initial
312 cases from the dataset were enrolled in the randomized trial. Thus,
the dataset specifically pertains to these 312 patients, for whom the
data is largely complete.

## Usage

``` r
data(pbc3)
```

## Format

A data frame with 1945 observations on the following 27 variables:

- `id`:

  patients identifier; in total there are 312 patients.

- `years`:

  number of years between registration and the earlier of death,
  transplantation, or study analysis time.

- `status`:

  a factor with levels `alive`, `transplanted` and `dead`.

- `drug`:

  a factor with levels `placebo` and `D-penicil`.

- `age`:

  at registration in years.

- `sex`:

  a factor with levels `male` and `female`.

- `year`:

  number of years between enrollment and this visit date, remaining
  values on the line of data refer to this visit.

- `ascites`:

  a factor with levels `No` and `Yes`.

- `hepatomegaly`:

  a factor with levels `No` and `Yes`.

- `spiders`:

  a factor with levels `No` and `Yes`.

- `edema`:

  a factor with levels `No edema` (i.e. no edema and no diuretic therapy
  for edema), `edema no diuretics` (i.e. edema present without
  diuretics, or edema resolved by diuretics), and
  `edema despite diuretics` (i.e. edema despite diuretic therapy).

- `serBilir`:

  serum bilirubin in mg/dl.

- `serChol`:

  serum cholesterol in mg/dl.

- `albumin`:

  albumin in mg/dl.

- `alkaline`:

  alkaline phosphatase in U/liter.

- `SGOT`:

  SGOT in U/ml.

- `platelets`:

  platelets per cubic ml/1000.

- `prothrombin`:

  prothrombin time in seconds.

- `histologic`:

  histologic stage of disease.

- `status2`:

  a numeric vector with the value 1 denoting if the patient was dead,
  and 0 if the patient was alive or transplanted.

- `status3`:

  a numeric vector with the value 1 denoting if the patient was dead or
  transplanted, and 0 if the patient was alive.

- `status4`:

  a numeric vector with the value 1 denoting if the patient was
  transplanted, and 0 if the patient was dead. Used for competing risks.

- `status5`:

  a numeric vector with the value 2 if the patient was transplanted, 1
  denoting if the patient was dead, and 0 if the patient was alive. Used
  for competing risks with censored.

- `Tyears1`:

  a numeric vector with a transformed value of time-to-event outcome.

- `Tyears2`:

  a numeric vector with a transformed value of time-to-event outcome.

- `Tyears3`:

  a numeric vector with a transformed value of time-to-event outcome.

- `Tyears4`:

  a numeric vector with a transformed value of time-to-event outcome.

## Source

[`pbc`](https://rdrr.io/pkg/survival/man/pbc.html).

## References

Fleming T, Harrington D. *Counting Processes and Survival Analysis*.
1991; New York: Wiley.

Therneau T, Grambsch P. *Modeling Survival Data: Extending the Cox
Model*. 2000; New York: Springer-Verlag.
