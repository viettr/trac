#' Two stage regression
#'
#' Perform a two stage fitting procedure similar to the one proposed by
#' Bates, S., & Tibshirani, R. (2019). Log‐ratio lasso:
#' Scalable, sparse estimation for log‐ratio models. Biometrics, 75(2), 613-624.
#' This function uses the selected components by trac or log-ratio and
#' fits a sparse model based on all possible log-ratios.
#'
#' 1. Fit trac or sparse-log contrast and extract the selected components.
#' 2. Fit a second-stage based on the selected components.
#'
#' @param Z n by p matrix containing log(X) (see \code{Z} from
#'     \code{\link{trac}})
#' @param y n vector (response) (see \code{y} from \code{\link{trac}})
#' @param A p by (t_size-1) binary matrix giving tree structure (t_size is the
#'   total number of nodes and the -1 is because we do not include the root).
#'   Only needed for trac based models. (see \code{A} from \code{\link{trac}}).
#'   If A = NULL, a sparse-log contrast model is assumed.
#' @param additional_covariates n by p' matrix containing additional covariates,
#'   (see \code{additional_covariates} from \code{\link{trac}}).
#'   Non-compositional components are currently not penalized by the lasso
#' @param betas pre-screened coefficients
#'   (see output \code{gamma} from \code{\link{trac}} or
#'   output \code{beta} from \code{\link{sparse_log_contrast}})
#' @param topk maximum number of pre-screened coefficients to consider.
#'   Default NULL.
#' @param nfolds number of folds
#' @param method string which estimation method to use should be "regr" or
#'   "classif"
#' @param criterion which criterion should be used to select the coefficients
#'   "1se" or "min"
#' @param alpha nudge the model to select on a higher or lower level of the
#'   tree. Only relevant for trac based models.
#' @param w_additional weight for the additional covaraites
#' @param folds predefined folds (see \code{\link{cv_trac}})
#' @param classo Should the solver c-lasso be used instead of glmnet? Usefull
#'   for smaller models
#'
#'
#' @return list with: log_ratios: betas for log ratios; index: dataframe with
#'   index of the pre selected coefficients and the log ratio name; A: taxonomic
#'   tree information; method: regression or classification; cv_glmnet:
#'   output of glmnet or c-lasso,
#'   useful for prediction; criterion: which criterion to
#'   be used to select lambda based on cv (cross validation)
#' @export

second_stage <- function (Z,
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
                          classo = FALSE) {
  criterion <- match.arg(criterion)
  method <- match.arg(method)
  stopifnot(alpha >= 0)
  if (!is.null(A)) {
    n_alphas_compositional <- ncol(A)
  }
  else {
    n_alphas_compositional <- ncol(Z)
  }
  pre_selected <- betas[1:n_alphas_compositional]
  if (is.null(names(pre_selected))) {
    names(pre_selected) <- paste0("gamma_", 1:length(pre_selected))
  }
  pre_selected_index <- abs(pre_selected) >= 1e-07
  if (length(pre_selected[pre_selected_index]) < 2) {
    stop(
      "Only two or less preselected variables are passed.\n
        Therefore using a second step doesn't make sense."
    )
  }
  if (!is.null(topk)) {
    if (!is.numeric(topk))
      stop("Topk has to be numeric")
    if (topk <= 2)
      stop("Minimum number of topk is 3")
    if (length(pre_selected[pre_selected_index]) > topk) {
      id <- 1:length(pre_selected)
      topk_df <-
        data.frame(
          id = id,
          pre_selected = abs(pre_selected),
          pre_selected_index = pre_selected_index
        )
      topk_df <-
        topk_df[order(topk_df$pre_selected, decreasing = TRUE), ]
      topk_df$pre_selected[(topk + 1):length(pre_selected)] <-
        0
      topk_df <-
        topk_df[order(topk_df$id, decreasing = FALSE), ]
      pre_selected <- topk_df$pre_selected
      names(pre_selected) <- rownames(topk_df)
      pre_selected_index <- abs(pre_selected) > 0
    }
  }
  if (!is.null(A)) {
    Z_A <- as.matrix(Z %*% A)
    subtree_n <- colSums(as.matrix(A))
    Z_A <- t((1 / subtree_n) * t(Z_A))
  }
  else {
    Z_A <- Z
  }
  index <-
    expand.grid(which(pre_selected_index), which(pre_selected_index))
  index <- index[index$Var2 > index$Var1, ]
  index <- as.data.frame(index)
  if (is.null(colnames(Z_A))) {
    colnames(Z_A) <- paste0("comp_", 1:length(pre_selected))
  }
  Z_A_names <- colnames(Z_A)
  n <- nrow(Z_A)
  colnames(index) <- c("index1", "index2")
  index$variable1 <- Z_A_names[index$index1]
  index$variable2 <- Z_A_names[index$index2]
  if (!is.null(A)) {
    penalty_weight_comp <-
      subtree_n[index$index1] * subtree_n[index$index2]
    penalty_weight_comp <- penalty_weight_comp ^ alpha
    penalty_weight_comp <- 1 / penalty_weight_comp
  }
  else {
    penalty_weight_comp <- rep(1, nrow(index))
  }
  index$log_name <- paste(index$variable1, index$variable2,
                          sep = "///")
  n_log_contrast <- nrow(index)
  expanded_z <- matrix(0, nrow = n, ncol = n_log_contrast)
  colnames(expanded_z) <- index$log_name
  for (i in 1:n_log_contrast) {
    expanded_z[, i] <- Z_A[, index[i, "index1"]] - Z_A[,
                                                       index[i, "index2"]]
  }
  n_comp <- ncol(expanded_z)
  if (!is.null(additional_covariates)) {
    categorical_list <-
      get_categorical_variables(additional_covariates)
    categorical <- categorical_list[["categorical"]]
    n_categorical <- categorical_list[["n_categorical"]]
    if (n_categorical > 0) {
      additional_covariates[, categorical] <-
        transform_categorical_variables(additional_covariates,
                                               categorical)
    }
    n_x <- ncol(additional_covariates)
    if (is.null(w_additional)) {
      penalty_factor <- c(penalty_weight_comp, rep(1, n_x))
    } else {
      penalty_factor <- c(penalty_weight_comp, rep(w_additional, n_x))
    }
    expanded_z <- cbind(expanded_z, additional_covariates)
    expanded_z <- as.matrix(expanded_z)
  }
  else {
    penalty_factor <- penalty_weight_comp
  }

  # transform list with folds to folds_ids containing a vector
  # indicating to which fold each observation belongs to
  if (classo == FALSE & !is.null(folds)) {
    fold_ids <- tibble::enframe(folds,
                                name = "name", value = "value") %>%
      # tibble with name as fold id and value as the position of the observation
      # therefore sort by value and extract only names
      tidyr::unnest(cols = c(value)) %>%
      dplyr::arrange("value") %>%
      dplyr::select("name") %>%
      c()
  }
  if (classo == TRUE) {
    fit <- classo_fitting(X = expanded_z, y = y, method = method)

    cvfit <- cv_classo_fitting(fit = fit, X = expanded_z, y = y,
                               folds = folds, nfolds = nfolds)
    expanded_z_normalized <-
      normalization_additional_covariates(
        additional_covariates = as.data.frame(as.matrix(expanded_z)),
        p_x = dim(expanded_z)[2], intercept = TRUE)



    log_ratios <- fit$beta[, cvfit$cv$ibest]

    if (criterion == "1se") {
      log_ratios <-  fit$beta[, cvfit$cv$i1se]
      beta0 <- fit$beta0[cvfit$cv$i1se]
    }
    if (criterion == "min") {
      log_ratios <- fit$beta[, cvfit$cv$ibest]
      beta0 <- fit$beta0[cvfit$cv$ibest]
    }
    return(list(
      log_ratios = log_ratios,
      beta0 = beta0,
      index = index,
      method = method,
      A = A,
      cv_glmnet = list(fit = fit, cvfit = cvfit),
      criterion = criterion,
      fit_method = "classo"))
  } else {

    if (method == "classif") {
      if (!is.factor(y)) {
        y_new <- as.factor(y)
      }
      else {
        y_new <- y
      }
    }
    else {
      y_new <- y
    }
    if (is.null(folds)) fold_ids <- NULL
    if (method == "regr") {
      if (is.null(fold_ids)) {
        cv_glmnet <- glmnet::cv.glmnet(
          x = expanded_z,
          y = y_new,
          alpha = 1,
          nfolds = nfolds,
          standardize = TRUE,
          type.measure = "mse",
          family = "gaussian",
          penalty.factor = penalty_factor
        )
      } else {
        cv_glmnet <- glmnet::cv.glmnet(
          x = expanded_z,
          y = y_new,
          alpha = 1,
          foldsid = fold_ids,
          standardize = TRUE,
          type.measure = "mse",
          family = "gaussian",
          penalty.factor = penalty_factor
        )
      }

    }
    else if (method == "classif") {
      if (is.null(fold_ids)) {
        cv_glmnet <- glmnet::cv.glmnet(
          x = expanded_z,
          y = y_new,
          alpha = 1,
          nfolds = nfolds,
          standardize = TRUE,
          family = "binomial",
          type.measure = "class",
          penalty.factor = penalty_factor
        )
      } else {
        cv_glmnet <- glmnet::cv.glmnet(
          x = expanded_z,
          y = y_new,
          alpha = 1,
          foldsid = fold_ids,
          standardize = TRUE,
          family = "binomial",
          type.measure = "class",
          penalty.factor = penalty_factor
        )
      }
    }
    if (criterion == "1se") {
      log_ratios <- as.matrix(glmnet::coef.glmnet(cv_glmnet,
                                                  s = "lambda.1se"))[, 1]
    }
    if (criterion == "min") {
      log_ratios <- as.matrix(glmnet::coef.glmnet(cv_glmnet,
                                                  s = "lambda.min"))[, 1]
    }

    return(list(
      log_ratios = log_ratios[2:length(log_ratios)],
      beta0 = log_ratios[1],
      index = index,
      A = A,
      method = method,
      cv_glmnet = cv_glmnet,
      criterion = criterion,
      fit_method = "glmnet"))
  }
}

#' Make predictions based on a second_stage fit
#'
#' @param new_Z a new data matrix (see \code{Z} from \code{\link{second_stage}})
#' @param new_additional_covariates a new data matrix
#'   (see \code{additional_covariates} from \code{\link{second_stage}})
#' @param fit output of the function \code{\link{second_stage}}
#' @param output string  either "raw", "probability" or "class" only relevant
#'   classification tasks
#' @return a vector of \code{nrow(new_Z) + nrow(new_additional_covariates)}
#'   predictions.
#' @importFrom stats predict
#' @export


predict_second_stage <- function(new_Z, new_additional_covariates = NULL, fit,
                                 output = c("raw", "probability", "class")) {
  # partial matching for input
  output <- match.arg(output)
  # create geometric mean
  if (!is.null(fit$A)) {
    A <- as.matrix(fit$A)
    Z_A <- as.matrix(new_Z %*% A)
    subtree_n <- colSums(A)
    Z_A <- t((1 / subtree_n) * t(Z_A))
  } else {
    Z_A <- new_Z
  }
  n <- nrow(Z_A)
  # create logratios
  index <- fit$index
  n_log_contrast <- nrow(index)
  expanded_z <- matrix(0, nrow = n, ncol = n_log_contrast)
  colnames(expanded_z) <- index$log_name
  for (i in 1:n_log_contrast) {
    expanded_z[, i] <- Z_A[, index[i, "index1"]] - Z_A[, index[i, "index2"]]
  }
  # add aditional covariates
  if (!is.null(new_additional_covariates)) {
    categorical_list <- get_categorical_variables(new_additional_covariates)
    categorical <- categorical_list[["categorical"]]
    n_categorical <- categorical_list[["n_categorical"]]
    if (n_categorical > 0) {
      new_additional_covariates[, categorical] <-
        transform_categorical_variables(new_additional_covariates, categorical)
    }
    expanded_z <- cbind(expanded_z, new_additional_covariates)
    expanded_z <- as.matrix(expanded_z)
  }

  if (is.null(fit$fit_method)) fit$fit_method <- "glmnet"

  if (fit$fit_method == "classo") {
    yhat <- predict_trac(list(fit$cv_glmnet$fit), new_Z = as.matrix(expanded_z))
    i_selected <- ifelse(fit$criterion == "1se", fit$cv_glmnet$cvfit$cv$i1se,
                         fit$cv_glmnet$cvfit$cv$ibest)
    yhat <- yhat[[1]][ , i_selected]
    return(yhat)
  } else {
    # specify how to select lambda
    if (fit$criterion == "1se") s <- "lambda.1se"
    if (fit$criterion == "min") s <- "lambda.min"
    # define different output possibilites for classification
    if (fit$method == "classif") {
      if (output == "raw") type <- "link"
      if (output == "probability") type <- "response"
      if (output == "class") type <- "class"
    } else if (fit$method == "regr") {
      type <- "link"
    }
    # make the prediction
    yhat <- predict(object = fit$cv_glmnet,
                    newx = expanded_z, s = s, type = type)
    # return prediction
    return(yhat)
  }
}
