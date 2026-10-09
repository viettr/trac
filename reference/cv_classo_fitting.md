# Perform cross validation for tuning parameter selection for sparse log contrast

This function is to be called after calling
[`sparse_log_contrast`](https://viettran.de/trac/reference/sparse_log_contrast.md).
It performs `nfold`-fold cross validation.

## Usage

``` r
cv_classo_fitting(
  fit,
  X,
  y,
  folds = NULL,
  nfolds = 5,
  summary_function = stats::median
)
```

## Arguments

- fit:

  output of
  [`sparse_log_contrast`](https://viettran.de/trac/reference/sparse_log_contrast.md)
  function.

- X, y:

  same arguments as passed to
  [`sparse_log_contrast`](https://viettran.de/trac/reference/sparse_log_contrast.md).
  C will be taken from fit object.

- folds:

  a partition of `1:nrow(Z)`.

- nfolds:

  number of folds for cross-validation

- summary_function:

  how to combine the errors calculated on each observation within a fold
  (e.g. mean or median)
