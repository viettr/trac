# Perform cross validation for tuning parameter selection

This function is to be called after calling
[`trac`](https://viettran.de/trac/reference/trac.md). It performs
`nfold`-fold cross validation. For classification the metric is
misclassification error.

## Usage

``` r
cv_trac(
  fit,
  Z,
  y,
  A,
  additional_covariates = NULL,
  folds = NULL,
  nfolds = 5,
  summary_function = stats::median
)
```

## Arguments

- fit:

  output of [`trac`](https://viettran.de/trac/reference/trac.md)
  function.

- Z, y, A, additional_covariates:

  same arguments as passed to
  [`trac`](https://viettran.de/trac/reference/trac.md)

- folds:

  a partition of `1:nrow(Z)`.

- nfolds:

  number of folds for cross-validation

- summary_function:

  how to combine the errors calculated on each observation within a fold
  (e.g. mean or median) (only for regression task)
