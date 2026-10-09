# Make predictions based on a sparse log contrast fit

Make predictions based on a sparse log contrast fit

## Usage

``` r
predict_sparse_log_contrast(
  fit,
  new_Z,
  new_additional_covariates = NULL,
  output = c("raw", "class")
)
```

## Arguments

- fit:

  output of the function
  [`sparse_log_contrast`](https://viettran.de/trac/reference/sparse_log_contrast.md)

- new_Z:

  a new data matrix (see `Z` from
  [`sparse_log_contrast`](https://viettran.de/trac/reference/sparse_log_contrast.md))

- new_additional_covariates:

  a new data matrix (see `additional_covariates` from
  [`sparse_log_contrast`](https://viettran.de/trac/reference/sparse_log_contrast.md))

- output:

  string either "raw" or "class" only relevant classification tasks

## Value

a vector of `nrow(new_Z)` predictions.
