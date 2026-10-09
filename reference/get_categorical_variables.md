# Check if additional variables are categorical

Check if the additional non-compositional covariates are categorical or
not based on the type of the column. If the column is a binary factor
then it is assumed to be a categorical variable. Useful to transform the
input.

## Usage

``` r
get_categorical_variables(X)
```

## Arguments

- X:

  see `additional_covariates` from
  [`trac`](https://viettran.de/trac/reference/trac.md)

## Value

list with vector indicating if categorical or not, number of categorical
variables
