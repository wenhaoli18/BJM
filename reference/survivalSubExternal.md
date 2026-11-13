# Survival sub-model that borrows a long-term tail from an external cohort

**Experimental.** Predicting beyond the last follow-up of the study
needs the baseline hazard there, which
[`survivalSub()`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md)
can only extrapolate (linearly in the cumulative hazard). When an
external cohort with longer follow-up is available – an earlier
generation of the same device, a registry, an older trial – its baseline
hazard can supply the shape of the tail. `survivalSubExternal()` fits
one Cox model to both cohorts with a generation indicator (1 for the
study, 0 for the external cohort): \$\$h(t \mid x, g) = h_0(t)
\exp(\beta^\top x + \gamma g),\$\$ so the study's hazard beyond its own
follow-up is the external one scaled by \\e^\gamma\\. Only the survival
sub-model uses the external cohort; the longitudinal sub-model is fit on
the study as usual, and the external cohort needs no biomarkers.

The key assumption is proportional hazards between the cohorts, also
beyond the study's follow-up, where it cannot be checked. Within the
overlap it is tested with
[`cox.zph`](https://rdrr.io/pkg/survival/man/cox.zph.html) (stored in
`$external$ph_test`, with a warning when \\p \< 0.05\\); the test has
little power with short study follow-up, so a non-significant result is
weak reassurance, and a sensitivity analysis on the tail is advisable.
The hazard ratio itself, estimated from the overlap only, can be noisy
when the study is short and has few events, and the proportional hazards
assumption carries that noise to the whole tail. `tail_hr` fixes it
instead (the other coefficients are then re-estimated given it), for a
sensitivity analysis over plausible values, e.g. `c(0.8, 1, 1.25)`, or
to use prior knowledge, e.g. that the generations are equivalent.

With `splice = TRUE` (the default) the external cohort only supplies the
tail: up to the study's last follow-up \\\tau\\ the baseline cumulative
hazard is the Breslow estimate from the study alone (at the coefficients
of the joint fit), and beyond \\\tau\\ it continues with the increments
of the joint fit's baseline, i.e. the external shape scaled by
\\e^\gamma\\. Within the study's follow-up predictions then rest on the
study's own events. With `splice = FALSE` the joint fit's baseline is
used throughout, so the external cohort (scaled by \\e^\gamma\\) also
shapes the study's hazard within its follow-up, which matters when the
cohorts are not proportional there.

Predictions are for the study population: when the prediction data have
no `generation` column it is taken as 1.

## Usage

``` r
survivalSubExternal(
  data_internal,
  data_external,
  form_marginal_surv,
  form_conditional_cr = NULL,
  generation = "generation",
  tail_hr = NULL,
  splice = TRUE
)
```

## Arguments

- data_internal, data_external:

  One row per subject, for the study and for the external cohort. Both
  need the variables of `form_marginal_surv` (and of
  `form_conditional_cr`, if given).

- form_marginal_surv, form_conditional_cr:

  As in
  [`survivalSub`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md),
  without the generation term, which is added.

- generation:

  Name of the generation indicator column to create (it must not already
  exist).

- tail_hr:

  `NULL` (default) to estimate the hazard ratio of the study against the
  external cohort, or a positive number to fix it.

- splice:

  If `TRUE` (default), use the study's own baseline up to its last
  follow-up and the external shape only beyond it (see Details). Not
  available for a stratified Cox model.

## Value

A `survivalSub.BJM` fit, usable wherever
[`survivalSub()`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md)'s
is (e.g.
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md),
`fitImputedBJM(survival_fit_all = )`), with an `external` element: the
generation variable, the cohort sizes and follow-up, the hazard ratio of
the study against the external cohort (`hazard_ratio`, the one used;
`hazard_ratio_estimated`, from the data; `tail_hr_fixed`) and the PH
test.

## Examples

``` r
# \donttest{
data(pbc3)
surv <- pbc3[!duplicated(pbc3$id), ]
# pretend the study had 4 years of follow-up and the rest of pbc3 is external
id <- as.numeric(as.character(surv$id))
study <- surv[id %% 2 == 0, ]
study$status3[study$years > 4] <- 0
study$years <- pmin(study$years, 4)
external <- surv[id %% 2 == 1, ]
fit <- survivalSubExternal(study, external, Surv(years, status3) ~ age + sex, NULL)
fit$external
#> $generation
#> [1] "generation"
#> 
#> $n_internal
#> [1] 156
#> 
#> $n_external
#> [1] 156
#> 
#> $max_followup_internal
#> [1] 4
#> 
#> $max_followup_external
#> [1] 14.30566
#> 
#> $hazard_ratio
#> [1] 0.6608038
#> 
#> $hazard_ratio_estimated
#> [1] 0.6608038
#> 
#> $tail_hr_fixed
#> [1] FALSE
#> 
#> $ph_test
#>     chisq        df         p 
#> 0.5637548 1.0000000 0.4527513 
#> 
#> $splice
#> [1] TRUE
#> 
#> $cum_basehaz
#>          hazard       time
#> 1   0.004250954  0.1396342
#> 2   0.008542223  0.1943927
#> 3   0.012859635  0.3614062
#> 4   0.017204817  0.3833096
#> 5   0.021617442  0.5092542
#> 6   0.026062732  0.5667506
#> 7   0.030540659  0.7556675
#> 8   0.035041742  0.8323294
#> 9   0.039577085  0.8788742
#> 10  0.044135043  1.0650531
#> 11  0.048746879  1.5031212
#> 12  0.053384210  1.5113350
#> 13  0.058071195  1.6345417
#> 14  0.062785369  1.8426240
#> 15  0.067534613  2.0178513
#> 16  0.072307761  2.0534443
#> 17  0.077114863  2.2560508
#> 18  0.081956141  2.3354507
#> 19  0.086837337  2.4367539
#> 20  0.091799793  2.4750849
#> 21  0.096806167  2.5873398
#> 22  0.101848459  2.6667397
#> 23  0.106916274  2.8392290
#> 24  0.112012732  2.9213668
#> 25  0.117143509  3.1540905
#> 26  0.122369900  3.2608696
#> 27  0.127638515  3.3320556
#> 28  0.132979213  3.3813383
#> 29  0.138373792  3.5510897
#> 30  0.143835322  3.7126273
#> 31  0.149348387  3.8632132
#> 32  0.154897293  3.9070200
#> 33  0.160472464  3.9261855
#> 34  0.166080082  3.9289234
#> 35  0.171715916  3.9535648
#> 36  0.177397016  3.9617786
#> 37  0.177397016  3.9891578
#> 38  0.177397016  4.0000000
#> 88  0.181779033  4.1205782
#> 89  0.181779033  4.1287920
#> 90  0.186220522  4.2054539
#> 91  0.186220522  4.2492608
#> 92  0.190773693  4.3724674
#> 93  0.195382145  4.5148396
#> 94  0.195382145  4.5531705
#> 95  0.200086060  4.6106670
#> 96  0.204888132  4.7667287
#> 97  0.204888132  4.9775490
#> 98  0.204888132  5.1609900
#> 99  0.209812247  5.1965831
#> 100 0.214780643  5.3307414
#> 101 0.214780643  5.3718103
#> 102 0.214780643  5.3909758
#> 103 0.214780643  5.4457343
#> 104 0.214780643  5.5306100
#> 105 0.220005629  5.6264374
#> 106 0.220005629  5.6675063
#> 107 0.225344509  5.6976235
#> 108 0.225344509  5.7359544
#> 109 0.225344509  5.7551199
#> 110 0.225344509  5.7633337
#> 111 0.225344509  5.9933195
#> 112 0.231003646  6.0891469
#> 113 0.231003646  6.1274778
#> 114 0.236783062  6.1356916
#> 115 0.236783062  6.2506845
#> 116 0.236783062  6.2534224
#> 117 0.242794894  6.2890154
#> 118 0.248891777  6.3027051
#> 119 0.248891777  6.3163947
#> 120 0.248891777  6.4286497
#> 121 0.255217218  6.4341255
#> 122 0.261609224  6.4532910
#> 123 0.268079175  6.5710218
#> 124 0.268079175  6.5847114
#> 125 0.268079175  6.6641113
#> 126 0.268079175  6.7133939
#> 127 0.268079175  6.7325594
#> 128 0.275170064  6.7572007
#> 129 0.275170064  6.8475523
#> 130 0.282480338  6.9543314
#> 131 0.282480338  7.1213449
#> 132 0.282480338  7.3677582
#> 133 0.282480338  7.4444201
#> 134 0.290272442  7.5813164
#> 135 0.298194292  7.6552404
#> 136 0.298194292  7.7373782
#> 137 0.298194292  7.7757091
#> 138 0.306780176  7.7948746
#> 139 0.306780176  7.8441573
#> 140 0.315683836  7.8907020
#> 141 0.315683836  7.9071296
#> 142 0.315683836  7.9126054
#> 143 0.315683836  8.1124740
#> 144 0.315683836  8.1891359
#> 145 0.315683836  8.2110393
#> 146 0.315683836  8.3068667
#> 147 0.325645581  8.4656664
#> 148 0.325645581  8.5560180
#> 149 0.335962592  8.6792246
#> 150 0.335962592  8.6874384
#> 151 0.335962592  8.8298105
#> 152 0.335962592  8.8763553
#> 153 0.335962592  8.8790932
#> 154 0.335962592  8.8845691
#> 155 0.335962592  8.9749206
#> 156 0.347996417  8.9858723
#> 157 0.347996417  9.0159895
#> 158 0.347996417  9.1063410
#> 159 0.360905413  9.1939547
#> 160 0.374139949  9.2596649
#> 161 0.374139949  9.4896506
#> 162 0.374139949  9.5936918
#> 163 0.374139949  9.6265469
#> 164 0.374139949  9.6703537
#> 165 0.374139949  9.6895192
#> 166 0.390247037  9.8127259
#> 167 0.390247037  9.9797393
#> 168 0.407764232 10.0125945
#> 169 0.407764232 10.0180703
#> 170 0.407764232 10.1823459
#> 171 0.428722393 10.3000767
#> 172 0.428722393 10.4561384
#> 173 0.451571363 10.5108969
#> 174 0.475322450 10.5492279
#> 175 0.475322450 10.6998138
#> 176 0.475322450 10.7436206
#> 177 0.475322450 11.1105027
#> 178 0.475322450 11.2008542
#> 179 0.475322450 11.6279706
#> 180 0.475322450 11.6361844
#> 181 0.475322450 11.9811631
#> 182 0.475322450 12.2905487
#> 183 0.475322450 12.5479137
#> 184 0.475322450 12.8518235
#> 185 0.475322450 12.9202716
#> 186 0.475322450 13.1310919
#> 187 0.475322450 13.2570365
#> 188 0.475322450 13.3035812
#> 189 0.475322450 13.3200088
#> 190 0.475322450 13.4185741
#> 191 0.475322450 13.4842843
#> 192 0.682366432 13.8922352
#> 193 0.682366432 14.3056620
#> 
# sensitivity to the hazard ratio between the cohorts
fits <- lapply(c(0.8, 1, 1.25), function(hr)
  survivalSubExternal(study, external, Surv(years, status3) ~ age + sex, NULL, tail_hr = hr))
# }
```
