# Check non-compositional inputs

Check if the additional non-compositional inputs have NAs and have the
same number of observations as the compositional inputs

## Usage

``` r
check_additional_covariates(
  additional_covariates,
  n,
  w_additional_covariates,
  p_x
)
```

## Arguments

- additional_covariates:

  new data matrix (see `additional_covariates` from
  [`trac`](https://viettran.de/trac/reference/trac.md))

- n:

  number of observations

- w_additional_covariates:

  weights for the estimation of the coefficients

- p_x:

  vector with number of additional non-compositional covariates

## Value

errors if the requirements are not met
