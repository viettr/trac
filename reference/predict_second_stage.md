# Make predictions based on a second_stage fit

Make predictions based on a second_stage fit

## Usage

``` r
predict_second_stage(
  new_Z,
  new_additional_covariates = NULL,
  fit,
  output = c("raw", "probability", "class")
)
```

## Arguments

- new_Z:

  a new data matrix (see `Z` from
  [`second_stage`](https://viettran.de/trac/reference/second_stage.md))

- new_additional_covariates:

  a new data matrix (see `additional_covariates` from
  [`second_stage`](https://viettran.de/trac/reference/second_stage.md))

- fit:

  output of the function
  [`second_stage`](https://viettran.de/trac/reference/second_stage.md)

- output:

  string either "raw", "probability" or "class" only relevant
  classification tasks and glmnet in second stage

## Value

a vector of `nrow(new_Z) + nrow(new_additional_covariates)` predictions.
