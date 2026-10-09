#' Make predictions based on a trac fit
#'
#' @param fit output of the function \code{\link{trac}}
#' @param new_Z a new data matrix (see \code{Z} from \code{\link{trac}})
#' @param new_additional_covariates a new data matrix
#'    (see \code{additional_covariates} from \code{\link{trac}})
#' @param output string  either "raw" or "class" only relevant for
#'   classification tasks
#' @return a vector of \code{nrow(new_Z)} predictions.
#' @export
predict_trac <- function(fit, new_Z, new_additional_covariates = NULL,
                         output = c("raw", "class")) {
  # fit: output of wag
  # new_Z: n_new by p matrix
  # new_additional_covariates: n_new by p' matrix



  # add to make make code backwards compatible, regression
  # and with intercept as default
  if (is.null(fit[[1]]$method)) {
    fit[[1]]$method <- "regr"
    fit[[1]]$intercept <- TRUE
  }
  # check which method to use
  output <- match.arg(output)
  classification <- fit[[1]]$method %in% c("classif",
                                           "classif_huber")
  intercept <- fit[[1]]$intercept
  
  # keep matrix if only one out of sample prediction and do not transform to
  # vector internally
  if (is.null(dim(new_Z)))                                                             # A14: one test row
    new_Z <- matrix(new_Z, nrow = 1, dimnames = list(NULL, names(new_Z)))
  # Transform additional covariates

  if (!is.null(new_additional_covariates)) {
    categorical_list <-
      get_categorical_variables(new_additional_covariates)
    categorical <- categorical_list[["categorical"]]
    n_categorical <- categorical_list[["n_categorical"]]
    if (n_categorical > 0) {
      new_additional_covariates[, categorical] <-
        transform_categorical_variables(new_additional_covariates,
                                        categorical)
    }
    new_Z <- cbind(as.matrix(new_Z), as.matrix(new_additional_covariates))
    new_Z <- as.matrix(new_Z)
  }

  # save output as list
  yhat <- list()
  for (iw in seq_along(fit)) {

    if (intercept) {
      yhat[[iw]] <- t(as.numeric(fit[[iw]]$beta0) +
        t(as.matrix(new_Z %*% fit[[iw]]$beta)))
    } else {
      yhat[[iw]] <- t(t(as.matrix(new_Z %*% fit[[iw]]$beta)))
    }
    # classification: transform output if not raw score
    if (classification) {
      if (output == "class") {
        yhat[[iw]] <- yhat[[iw]] >= 0
        yhat[[iw]] <- yhat[[iw]] * 2 - 1
      }
    }

    rownames(yhat[[iw]]) <- rownames(new_Z)
  }
  return(yhat)
}
