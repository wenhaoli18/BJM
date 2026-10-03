# Factor one patient's joint longitudinal covariance

The stacked (across every biomarker) observation vector of one patient
has covariance `Sigma_all = A %*% D %*% t(A) + diag(r_diag)`: `A` is the
`N x r` random-effects design
([`random_effects_design()`](https://wenhaoli18.github.io/BJM/reference/random_effects_design.md)),
`D` the `r x r` random-effects covariance (`Sigma_fit`) and `r_diag`
each observation's residual variance. This factors it once so that
[`cov_logdens()`](https://wenhaoli18.github.io/BJM/reference/cov_logdens.md)
can evaluate the Gaussian log-density at any number of residual vectors.

With `method = "woodbury"`, `Sigma_all` is never formed. Writing
`D = G %*% t(G)` and `B = A %*% G`, the matrix determinant lemma and the
Woodbury identity give \$\$\log\|\Sigma\| = \sum \log r_j +
\log\|K\|,\quad e'\Sigma^{-1}e = e'R^{-1}e - u'K^{-1}u,\$\$ with
`K = I + t(B) %*% R^{-1} %*% B` (`r x r`) and
`u = t(B) %*% R^{-1} %*% e`. Factoring costs `O(N r^2 + r^3)` instead of
`O(N^3)`, each density `O(N r)` instead of `O(N^2)`, and memory `O(N r)`
instead of `O(N^2)`. `D^{-1}` is never needed: `G` is `t(chol(D))`, or,
when `D` is only positive semi-definite (a random-effect variance
collapsing to 0 during EM), its eigenvectors scaled by the square roots
of the non-negative eigenvalues, which still reproduces `D` exactly.

With `method = "dense"`, `Sigma_all` is formed and Cholesky factored, as
[`mvtnorm::dmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/Mvnorm.html)
does. `"auto"` picks `"woodbury"` only when every residual variance is
positive and `N >= 3 r`: in timings the two break even around
`N = 2.5 r` (the Woodbury steps carry larger constants), and dense is
faster below. Both give the same density up to floating-point rounding.

## Usage

``` r
cov_factor(A, D, r_diag, method = c("auto", "woodbury", "dense"))
```

## Arguments

- A:

  Random-effects design matrix, `N x r`.

- D:

  Random-effects covariance matrix, `r x r`.

- r_diag:

  Residual variances, length `N`.

- method:

  `"auto"`, `"woodbury"` or `"dense"`.

## Value

An object for
[`cov_logdens()`](https://wenhaoli18.github.io/BJM/reference/cov_logdens.md).
