# Perform sparse log-contrast regression

Solves the constrained lasso problem using the CLASSO module in Python.
The optimization problem are

## Usage

``` r
sparse_log_contrast(
  Z,
  y,
  additional_covariates = NULL,
  C = NULL,
  fraclist = NULL,
  nlam = 20,
  min_frac = 1e-04,
  method = c("regr", "classif", "classif_huber"),
  w_additional_covariates = NULL,
  intercept = TRUE,
  rho = 0,
  limit_active = TRUE
)
```

## Arguments

- Z:

  n by p matrix containing log(X)

- y:

  n vector (response)

- additional_covariates:

  n by p' matrix containing additional covariates / features

- C:

  m by p matrix. Default is a row vector of ones.

- fraclist:

  (optional) vector of tuning parameter multipliers. Should be in (0,
  1\].

- nlam:

  number of tuning parameters (ignored if fraclist non-NULL)

- min_frac:

  smallest value of tuning parameter multiplier (ignored if fraclist
  non-NULL)

- method:

  string which estimation method to use should be in ("regr", "classif",
  "classif_huber")

- w_additional_covariates:

  vector of positive weights of length ncol(additional_covariates)
  (default: all equal to 1).

- intercept:

  only works for classification! Should the intercept be fitted. Default
  is TRUE, set to FALSE if the intercept should not be included

- rho:

  value for huberized classification loss. Default = 0.0.

- limit_active:

  (default = TRUE) should the solver stop optimizing if there are more
  non-zero features than degrees of freedom (in this case the number of
  observations)

## Details

Regression: minimize_beta, beta0 1/(2n) \|\| y - beta0 1_n - Zagg_clr
beta \|\|^2 + lamda_max \* frac \|\| beta \|\|*1 subject to C beta = 0
Classification: minimize*beta, beta0 max(1 - y_i(beta0 + Z_clr_i \*
beta), 0)^2 + lambda_max \* frac \|\| W \* beta \|\|\_1 subject to C
beta = 0

Default is C = 1_p^T, but C can be a general matrix.

Observe that the tuning parameter is specified through "frac", the
fraction of lamda_max (which is the smallest value for which beta is
nonzero).
