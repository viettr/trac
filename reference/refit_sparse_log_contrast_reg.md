# Refit log-contrast regression to sparsity constraint

Given output of
[`sparse_log_contrast`](https://viettran.de/trac/reference/sparse_log_contrast.md),
solves the regression problem with compositional constraint on the
features selected by sparse_log_contrast. In contrast to
[`refit_sparse_log_contrast`](https://viettran.de/trac/reference/refit_sparse_log_contrast.md)
this function only does the refit for one model specified by i_selected
or component_selected. This i_selected usually comes from the i1se from
the cross-validation output.

## Usage

``` r
refit_sparse_log_contrast_reg(
  fit,
  i_selected = NULL,
  Z,
  y,
  additional_covariates = NULL,
  tol = 1e-05,
  component_selected = NULL
)
```

## Arguments

- fit:

  output of sparse_log_contrast

- i_selected:

  indicator which lambda is selected based on for example the
  cross-validation procedure

- Z, y, additional_covariates:

  same arguments as passed to
  [`sparse_log_contrast`](https://viettran.de/trac/reference/sparse_log_contrast.md)

- tol:

  tolerance for deciding whether a beta value is zero

- component_selected:

  vector with indices which component to include
