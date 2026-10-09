# Make predictions based on a trac fit

Make predictions based on a trac fit

## Usage

``` r
predict_trac(
  fit,
  new_Z,
  new_additional_covariates = NULL,
  output = c("raw", "class")
)
```

## Arguments

- fit:

  output of the function
  [`trac`](https://viettran.de/trac/reference/trac.md)

- new_Z:

  a new data matrix (see `Z` from
  [`trac`](https://viettran.de/trac/reference/trac.md))

- new_additional_covariates:

  a new data matrix (see `additional_covariates` from
  [`trac`](https://viettran.de/trac/reference/trac.md))

- output:

  string either "raw" or "class" only relevant for classification tasks

## Value

a vector of `nrow(new_Z)` predictions.
