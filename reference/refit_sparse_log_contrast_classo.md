# Refit log-contrast regression/classification to sparsity constraint

Given output of
[`sparse_log_contrast`](https://viettran.de/trac/reference/sparse_log_contrast.md),
solves the regression problem with compositional constraint on the
features selected by sparse_log_contrast. In contrast to
[`refit_sparse_log_contrast`](https://viettran.de/trac/reference/refit_sparse_log_contrast.md)
this function refit the models based on c-lasso. Be aware that this
procedure does not guard against sign switches meaning that the signs of
the coefficients can change during the refitting procedure.

## Usage

``` r
refit_sparse_log_contrast_classo(
  fit,
  Z,
  y,
  additional_covariates = NULL,
  tol = 1e-05
)
```

## Arguments

- fit:

  output of sparse_log_contrast

- Z, y, additional_covariates:

  same arguments as passed to
  [`sparse_log_contrast`](https://viettran.de/trac/reference/sparse_log_contrast.md)

- tol:

  tolerance for deciding whether a beta value is zero
