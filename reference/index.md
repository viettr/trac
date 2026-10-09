# Package index

## Preprocessing

Functions for preparing your data for `trac`.

- [`aggregate_to_level()`](https://viettran.de/trac/reference/aggregate_to_level.md)
  : Aggregate data to a fixed level
- [`tax_table_to_phylo()`](https://viettran.de/trac/reference/tax_table_to_phylo.md)
  : Convert a tax table to a phylo object
- [`phylo_to_A()`](https://viettran.de/trac/reference/phylo_to_A.md) :
  Convert from phylo to the A matrix used in trac. Note this is similar
  to the A used in rare, but with the column of all ones (for the root)
  removed

## Main trac functionality

Functions for fitting `trac`, choosing tuning parameter, and making
predictions

- [`trac()`](https://viettran.de/trac/reference/trac.md) : Perform
  tree-based aggregation
- [`cv_trac()`](https://viettran.de/trac/reference/cv_trac.md) : Perform
  cross validation for tuning parameter selection
- [`predict_trac()`](https://viettran.de/trac/reference/predict_trac.md)
  : Make predictions based on a trac fit

## Plotting functions

Functions for vizualization of the fitting process

- [`plot_trac_path()`](https://viettran.de/trac/reference/plot_trac_path.md)
  : Plot trac coefficient path
- [`plot_cv_trac()`](https://viettran.de/trac/reference/plot_cv_trac.md)
  : Make a plot of the output of cv_trac

## Sparse log contrast functions

- [`sparse_log_contrast()`](https://viettran.de/trac/reference/sparse_log_contrast.md)
  : Perform sparse log-contrast regression
- [`cv_sparse_log_contrast()`](https://viettran.de/trac/reference/cv_sparse_log_contrast.md)
  : Perform cross validation for tuning parameter selection for sparse
  log contrast
- [`predict_sparse_log_contrast()`](https://viettran.de/trac/reference/predict_sparse_log_contrast.md)
  : Make predictions based on a sparse log contrast fit

- [`plot_trac_path()`](https://viettran.de/trac/reference/plot_trac_path.md)
  : Plot trac coefficient path
- [`plot_cv_trac()`](https://viettran.de/trac/reference/plot_cv_trac.md)
  : Make a plot of the output of cv_trac

## Example data sets

Gut microbiome data sets used in vignette

- [`sCD14`](https://viettran.de/trac/reference/sCD14.md) : sCD14 data
- [`malawi`](https://viettran.de/trac/reference/malawi.md) : malawi vs
  venezuela, adults only

## Refitting functions

Functions for fitting unregularized versions subject to the sparsity
constraints learned by trac or sparse log-contrast

- [`refit_trac()`](https://viettran.de/trac/reference/refit_trac.md) :
  Refit subject to sparsity constraints
- [`refit_sparse_log_contrast()`](https://viettran.de/trac/reference/refit_sparse_log_contrast.md)
  : Refit subject to sparsity constraints
- [`refit_sparse_log_contrast_classif()`](https://viettran.de/trac/reference/refit_sparse_log_contrast_classif.md)
  : Refit log-contrast for classification to sparsity constraint

## Extention: Two-stage procedure

Functions for fitting a second stage

- [`second_stage()`](https://viettran.de/trac/reference/second_stage.md)
  : Two stage regression
- [`predict_second_stage()`](https://viettran.de/trac/reference/predict_second_stage.md)
  : Make predictions based on a second_stage fit

## Additional topics

- [`classo_fitting()`](https://viettran.de/trac/reference/classo_fitting.md)
  : Perform sparse log-contrast regression
- [`cps_moisture`](https://viettran.de/trac/reference/cps_moisture.md) :
  Central Park soil moisture data
- [`cv_classo_fitting()`](https://viettran.de/trac/reference/cv_classo_fitting.md)
  : Perform cross validation for tuning parameter selection for sparse
  log contrast
- [`refit_sparse_log_contrast_classo()`](https://viettran.de/trac/reference/refit_sparse_log_contrast_classo.md)
  : Refit log-contrast regression/classification to sparsity constraint
- [`refit_sparse_log_contrast_reg()`](https://viettran.de/trac/reference/refit_sparse_log_contrast_reg.md)
  : Refit log-contrast regression to sparsity constraint
