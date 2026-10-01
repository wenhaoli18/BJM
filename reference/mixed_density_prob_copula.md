# Joint density x probability for a mixed continuous/ordinal stacked vector

Given one patient's full stacked mean vector and covariance matrix
(across all jointly-fit biomarkers' repeated measurements), computes the
exact Gaussian density at the observed values for the continuous-marker
rows, times the Gaussian-copula probability that the unobserved
ordinal-marker latent scores fall in the box implied by their observed
categories' cumulative-link thresholds – conditional on the continuous
rows' observed values (via the standard multivariate-normal
conditioning/Schur-complement formula). This is the quantity
[`conditionalYT()`](https://wenhaoli18.github.io/BJM/reference/conditionalYT.md)/[`conditionalYDT()`](https://wenhaoli18.github.io/BJM/reference/conditionalYDT.md)
compute directly as a single Gaussian density when every marker is
continuous (see the quadratic-form comment in `R/conditionalYT.R`);
setting `oo_idx = integer(0)` here reduces to exactly that computation
(the full `(2*pi)^{-n/2}` normalizing constant is kept, not dropped,
even though it would cancel in
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)'s
own ratio –
[`dynamicPredictionBio()`](https://wenhaoli18.github.io/BJM/reference/dynamicPredictionBio.md)
instead divides this function's output by
[`mvtnorm::dmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/Mvnorm.html)-normalized
values from
[`conditionalYTBioCopula()`](https://wenhaoli18.github.io/BJM/reference/conditionalYTBioCopula.md)/[`conditionalYDTBioCopula()`](https://wenhaoli18.github.io/BJM/reference/conditionalYDTBioCopula.md),
so the constant must be correct in absolute terms, not just consistent
within one call, for `Y_density` to remain a genuine, correctly-scaled
conditional density – see `R/conditionalYTBio.R`'s analogous
all-continuous computation, which keeps its own constant for the same
reason). Returns the **log** of that density/probability, so that it
neither overflows nor underflows (see
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
for how it is exponentiated); `-Inf` if the ordinal box probability is
0.

## Usage

``` r
mixed_density_prob_copula(
  Sigma_all,
  mu_full,
  y_all_vec,
  alpha_lower_vec,
  alpha_upper_vec,
  cc_idx,
  oo_idx
)
```

## Arguments

- Sigma_all:

  Full stacked covariance matrix for this patient (as returned by
  [`build_conditional_design_copula()`](https://wenhaoli18.github.io/BJM/reference/build_conditional_design_copula.md)).

- mu_full:

  Full stacked mean vector, same row order as `Sigma_all`.

- y_all_vec:

  Full stacked observed-value vector; entries at `oo_idx` are ignored.

- alpha_lower_vec, alpha_upper_vec:

  Full stacked threshold-bound vectors; entries at `cc_idx` are ignored.

- cc_idx, oo_idx:

  Integer index vectors (into the stacked row order) of the continuous-
  and ordinal-marker rows respectively.
