# Refit subject to sparsity constraints

Given output of [`trac`](https://viettran.de/trac/reference/trac.md),
solves the least squares problem with compositional constraint on the
features selected by trac.

## Usage

``` r
refit_trac(fit, Z, y, A, tol = 1e-05)
```

## Arguments

- fit:

  output of trac

- Z, y, A:

  same arguments as passed to
  [`trac`](https://viettran.de/trac/reference/trac.md)

- tol:

  tolerance for deciding whether a gamma value is zero

## Details

minimize_beta, beta0, gamma 1/(2n) \|\| y - beta0 1_n - Z beta \|\|^2
subject to gamma_trac nonselected = 0, beta = A gamma, 1_p^T beta = 0
