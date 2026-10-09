# Refit subject to sparsity constraints

Given output of
[`sparse_log_contrast`](https://viettran.de/trac/reference/sparse_log_contrast.md),
solves the least squares problem with compositional constraint on the
features selected by sparse_log_contrast.

## Usage

``` r
refit_sparse_log_contrast(fit, Z, y, tol = 1e-05)
```

## Arguments

- fit:

  output of sparse_log_contrast

- Z, y:

  same arguments as passed to
  [`sparse_log_contrast`](https://viettran.de/trac/reference/sparse_log_contrast.md)

- tol:

  tolerance for deciding whether a beta value is zero

## Details

minimize_beta, beta0 1/(2n) \|\| y - beta0 1_n - Z beta \|\|^2 subject
to beta_nonselected = 0, 1_p^T beta = 0
