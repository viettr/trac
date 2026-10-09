# Two stage regression

Perform a two stage fitting procedure similar to the one proposed by
Bates, S., & Tibshirani, R. (2019). Log‐ratio lasso: Scalable, sparse
estimation for log‐ratio models. Biometrics, 75(2), 613-624. This
function uses the selected components by trac or log-ratio and fits a
sparse model based on all possible log-ratios.

## Usage

``` r
second_stage(
  Z,
  A = NULL,
  y,
  additional_covariates = NULL,
  betas,
  topk = NULL,
  nfolds = 5,
  method = c("regr", "classif"),
  criterion = c("1se", "min"),
  alpha = 0,
  w_additional = NULL,
  folds = NULL,
  classo = FALSE,
  disjoint_filter = FALSE
)
```

## Arguments

- Z:

  n by p matrix containing log(X) (see `Z` from
  [`trac`](https://viettran.de/trac/reference/trac.md))

- A:

  p by (t_size-1) binary matrix giving tree structure (t_size is the
  total number of nodes and the -1 is because we do not include the
  root). Only needed for trac based models. (see `A` from
  [`trac`](https://viettran.de/trac/reference/trac.md)). If A = NULL, a
  sparse-log contrast model is assumed.

- y:

  n vector (response) (see `y` from
  [`trac`](https://viettran.de/trac/reference/trac.md))

- additional_covariates:

  n by p' matrix containing additional covariates, (see
  `additional_covariates` from
  [`trac`](https://viettran.de/trac/reference/trac.md)).
  Non-compositional components are currently not penalized by the lasso

- betas:

  pre-screened coefficients (see output `gamma` from
  [`trac`](https://viettran.de/trac/reference/trac.md) or output `beta`
  from
  [`sparse_log_contrast`](https://viettran.de/trac/reference/sparse_log_contrast.md))

- topk:

  maximum number of pre-screened coefficients to consider. Default NULL.

- nfolds:

  number of folds

- method:

  string which estimation method to use should be "regr" or "classif"

- criterion:

  which criterion should be used to select the coefficients "1se" or
  "min"

- alpha:

  nudge the model to select on a higher or lower level of the tree. Only
  relevant for trac based models.

- w_additional:

  weight for the additional covariates, does not work for c-lasso as
  stage II

- folds:

  predefined folds (see
  [`cv_trac`](https://viettran.de/trac/reference/cv_trac.md))

- classo:

  Should the solver c-lasso be used instead of glmnet? Usefull for
  smaller models

- disjoint_filter:

  If the preselected coefficients live on a tree we can enforce that no
  log-ratios within the same branch are selected

## Value

list with: log_ratios: betas for log ratios; index: dataframe with index
of the pre selected coefficients and the log ratio name; A: taxonomic
tree information; method: regression or classification; cv_glmnet:
output of glmnet or c-lasso, useful for prediction; criterion: which
criterion to be used to select lambda based on cv (cross validation).
Returns NULL with a warning if fewer than 2 variables survive
pre-screening or filtering, or if only one log-ratio and no additional
covariates remain (e.g. two parts, or a single disjoint pair after
`disjoint_filter`).

## Details

1.  Fit trac or sparse-log contrast and extract the selected components.

2.  Fit a second-stage based on the selected components.
