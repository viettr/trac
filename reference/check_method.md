# Check method input and other hyperparameter

Check if the method input and other hyperparameter for classification
are correctly specified

## Usage

``` r
check_method(method, y, rho = 0)
```

## Arguments

- method:

  The method (see `method` from
  [`trac`](https://viettran.de/trac/reference/trac.md))

- y:

  The outcome (see `y` from
  [`trac`](https://viettran.de/trac/reference/trac.md))

- rho:

  The hyperparameter for huberized loss for classification (see
  `rho`from [`trac`](https://viettran.de/trac/reference/trac.md))

## Value

list with vector indicating which method is used and the outcome
