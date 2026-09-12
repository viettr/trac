# Transform non-compositional input --------------------------------------------

#' Check if additional variables are categorical
#'
#' Check if the additional non-compositional covariates are categorical or not
#' based on the type of the column. If the column is a binary factor then it is
#' assumed to be a categorical variable. Useful to transform the input.
#'
#' @keywords internal
#' @param X see \code{additional_covariates} from \code{\link{trac}}
#' @return list with vector indicating if categorical or not, number of
#'    categorical variables
#'

get_categorical_variables <- function(X) {
# fct to get categorical covariates
  if (!is.data.frame(X)) X <- as.data.frame(X)
  categorical <- vapply(X, is.factor, logical(1))
  n_categorical <- sum(categorical)
  list(categorical = categorical, n_categorical = n_categorical)
}

#' Transform non-compositional categorical input
#'
#' Transform factors to numerical values of 0 and 1
#' @keywords internal
#' @param X see \code{additional_covariates} from \code{\link{trac}}
#' @param categorical vector indicating if categorical or not
#' @return data frame with transformed categorical variables from factor to
#'    numerical
#'

transform_categorical_variables <- function(X, categorical) {
  X <- as.data.frame(X)
  as.data.frame(lapply(X[categorical], function(x) as.numeric(x) - 1))
}

#' Check non-compositional inputs
#'
#' Check if the additional non-compositional inputs have NAs and have the same
#' number of observations as the compositional inputs
#' @keywords internal
#' @param additional_covariates new data matrix
#'    (see \code{additional_covariates} from \code{\link{trac}})
#' @param n number of observations
#' @param w_additional_covariates weights for the estimation of the coefficients
#' @param p_x vector with number of additional non-compositional covariates
#' @return errors if the requirements are not met
#' @export
#'


check_additional_covariates <-
  function(additional_covariates, n, w_additional_covariates, p_x) {
    # basic check input additional covariates
    stopifnot(nrow(additional_covariates) == n)
    # No missing data allowed
    if (any(is.na(additional_covariates))) {
      stop(paste(
        "missing data is currently not supported.",
        "There seems to be missing values in the non compositional covariates ",
        "(additional_covariates).",
        sep = " "
      ))
    }
    # Weight vector needs to be the same length as the number of covariates
    if (!is.vector(w_additional_covariates)) {
      stop("w_additional_covariates must be a matrix or a vector.")
    }
    if (length(w_additional_covariates) != p_x) {
      stop("w_additional_covariates must be
           of length ncol(additional_covariates)")
    }
}


# Check Input for classification -----------------------------------------------

#' Check method input and other hyperparameter
#'
#' Check if the method input and other hyperparameter for classification
#' are correctly specified
#'
#' @keywords internal
#' @param method The method (see \code{method}
#'     from \code{\link{trac}})
#' @param y The outcome (see \code{y}
#'     from \code{\link{trac}})
#' @param rho The hyperparameter for huberized loss for classification
#'     (see \code{rho}from \code{\link{trac}})
#' @return list with vector indicating which method is used and the outcome
#' @export
#'

check_method <- function(method, y, rho = 0.0) {
  # check the inputs for classification tasks
  supported_methods <- c("regr", "classif", "classif_huber")
  if (!(method %in% supported_methods)) {
    stop(paste("trac currently supports the following methods: ",
               paste(supported_methods, collapse = ", "),
               sep = " "
    ))
  }
  # If the outcome is binary warn the user
  if (length(unique(y)) == 2 & method == "regr") {
    warning("this looks like a classification task, check the method argument")
  }
  # create a dummy variable with information about weather it is a
  # classification task or not for later
  if (method %in% c("classif", "classif_huber")) {
    classification <- TRUE
  } else {
    classification <- FALSE
  }
  # The response variable should be binary for classification
  if (length(unique(y)) != 2 & classification) {
    stop("the response variable should be binary for the classification task")
  }
  # Transform the output for classification to 1, -1 if not done yet
  if (!all(unique(y) %in% c(-1, 1)) & classification) {
    stop("for classification y must be numeric with values -1 and 1")
  }
  # If not a classification task --> intercept is TRUE for
  # prediction
  # only accept rho for huberized loss smaller 1
  stopifnot(rho < 1)
  # return classification indicator and transformed y
  list(classification = classification,
       y = y)
}
