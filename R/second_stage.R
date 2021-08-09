#' Two stage regression
#'
#' @param Z n by p matrix containing log(X)
#' @param y n vector (response)
#' @param A p by (t_size-1) binary matrix giving tree structure (t_size is the
#'   total number of nodes and the -1 is because we do not include the root)
#' @param X n by p' matrix containing additional covariates
#' @param betas pre screened coefficients
#' @param topk maximum number of pre screened coefficients to consider.
#'   Default NULL.
#' @param nfolds number of folds
#' @param method string which estimation method to use should be "regr" or
#'   "classif"
#' @param criterion which criterion should be used to select the coefficients
#'   "1se" or "min"
#' @param alpha
#' @return list with: log_ratios: betas for log ratios; index: dataframe with
#'   index of the pre selected coefficients and the log ratio name; A: taxonomic
#'   tree information; method: regression or classification; cv_glmnet:
#'   output of glmnet, useful for prediction; criterion: which criterion to
#'   be used to select lambda based on cv (cross validation)
#' @export

second_stage <- function(Z, A, y, X = NULL, betas, topk = NULL, nfolds = 5,
                         method = c("regr", "classif"),
                         criterion = c("1se", "min"),
                         alpha = 0) {

  criterion <- match.arg(criterion)
  method <- match.arg(method)
  stopifnot(alpha >= 0)

  n_alphas_compositional <- ncol(A)
  pre_selected <- betas[1:n_alphas_compositional]
  if (is.null(names(pre_selected))) {
    names(pre_selected) <- paste0("gamma_", 1:length(pre_selected))
  }

  pre_selected_index <- abs(pre_selected) >= 1e-3

  if(length(pre_selected[pre_selected_index]) <= 2) {
    stop("Only two or less preselected variables are passed.
       Therefore using a second step doesn't make sense.")
  }

  if (!is.null(topk)) {
    if (!is.numeric(topk)) stop("Topk has to be numeric")
    if (topk <= 2) stop ("Minimum number of topk is 3")
    if (length(pre_selected[pre_selected_index]) > topk) {
      # create df in order to get topk values and replace the rest with 0
      # reorder the df into original order
      id <- 1:length(pre_selected)
      topk_df <- data.frame(id = id, pre_selected = abs(pre_selected),
                            pre_selected_index = pre_selected_index)
      topk_df <- topk_df[order(topk_df$pre_selected, decreasing = TRUE) ,]
      topk_df$pre_selected[(topk + 1):length(pre_selected)] <- 0
      topk_df <- topk_df[order(topk_df$id, decreasing = FALSE) ,]
      pre_selected <- topk_df$pre_selected
      names(pre_selected) <- rownames(topk_df)
      pre_selected_index <- abs(pre_selected) > 0
    }
  }


  # Change to geometric mean --> see formula 2

  Z_A <- as.matrix(Z %*% A)
  subtree_n <- colSums(as.matrix(A))
  Z_A <- t((1 / subtree_n) * t(Z_A))
#  Z_A <- scale(Z_A, center = TRUE, scale = FALSE)
  index <- expand.grid(which(pre_selected_index), which(pre_selected_index))
  index <- index[index$Var2 > index$Var1, ]
  # choose(length(pre_selected[pre_selected_index]), 2)
  index <- as.data.frame(index)
  Z_A_names <- colnames(Z_A)
  n <- nrow(Z_A)
  colnames(index) <- c("index1", "index2")

  index$variable1 <- Z_A_names[index$index1]
  index$variable2 <- Z_A_names[index$index2]

  penalty_weight_comp <- subtree_n[index$index1] * subtree_n[index$index2]
  penalty_weight_comp <- penalty_weight_comp^alpha
  penalty_weight_comp <- 1 / penalty_weight_comp


  index$log_name <- paste(index$variable1,
                          index$variable2, sep = "///")
  n_log_contrast <- nrow(index)


  expanded_z <- matrix(0, nrow = n, ncol = n_log_contrast)
  colnames(expanded_z) <- index$log_name

  for (i in 1:n_log_contrast) {
    expanded_z[, i] <- Z_A[, index[i, "index1"]] - Z_A[, index[i, "index2"]]
  }
  n_comp <- ncol(expanded_z)


  if (!is.null(X)) {
    categorical_list <- get_categorical_variables(X)
    categorical <- categorical_list[["categorical"]]
    n_categorical <- categorical_list[["n_categorical"]]
    if (n_categorical > 0) {
      X[, categorical] <- sapply(X[, categorical], as.numeric)
      X[, categorical] <- sapply(X[, categorical], function(x) x - 1)
    }
    n_x <- ncol(X)
    penalty_factor <- c(penalty_weight_comp, rep(0, n_x))
    expanded_z <- cbind(expanded_z, X)
    expanded_z <- as.matrix(expanded_z)
  } else {
    penalty_factor <- penalty_weight_comp
  }


  if (method == "classif") {
    if (!is.factor(y)) {
      y_new <- as.factor(y)
    } else {
      y_new <- y
    }
  }

  if (method == "regr") {
    cv_glmnet <- glmnet::cv.glmnet(x = expanded_z, y = y_new, alpha = 1,
                                   nfolds = nfolds, standardize = TRUE,
                                   family = "gaussian", penalty.factor =
                                     penalty_factor)
  } else if (method == "classif") {
    cv_glmnet <- glmnet::cv.glmnet(x = expanded_z, y = y_new, alpha = 1,
                                   nfolds = nfolds, standardize = TRUE,
                                   family = "binomial", type.measure = "class",
                                   penalty.factor = penalty_factor)
  }

  if (criterion == "1se") {
    log_ratios <- as.matrix(
      glmnet::coef.glmnet(cv_glmnet, s = "lambda.1se"))[, 1]
  }
  if (criterion == "min") {
    log_ratios <- as.matrix(
      glmnet::coef.glmnet(cv_glmnet, s = "lambda.min"))[, 1]
  }



  list(log_ratios = log_ratios[2:length(log_ratios)],
       beta0 = log_ratios[1],
       index = index,
       A = A,
       method = method,
       cv_glmnet = cv_glmnet,
       criterion = criterion)
}

#' Make predictions based on a second_stage fit
#'
#' @param new_Z a new data matrix (see \code{Z} from \code{\link{second_stage}})
#' @param new_X a new data matrix (see \code{X} from \code{\link{second_stage}})
#' @param fit output of the function \code{\link{second_stage}}
#' @param output string  either "raw", "probability" or "class" only relevant
#'   classification tasks
#' @return a vector of \code{nrow(new_Z) + nrow(new_X)} predictions.
#' @export


predict_second_stage <- function(new_Z, new_X = NULL, fit,
                                 output = c("raw", "probability", "class")) {
  output <- match.arg(output)

  A <- as.matrix(fit$A)
  Z_A <- as.matrix(new_Z %*% A)
  subtree_n <- colSums(A)
  Z_A <- t((1 / subtree_n) * t(Z_A))
  n <- nrow(Z_A)

  index <- fit$index
  n_log_contrast <- nrow(index)

  expanded_z <- matrix(0, nrow = n, ncol = n_log_contrast)
  colnames(expanded_z) <- index$log_name

  for (i in 1:n_log_contrast) {
    expanded_z[, i] <- Z_A[, index[i, "index1"]] - Z_A[, index[i, "index2"]]
  }
  if (!is.null(new_X)) {
    categorical_list <- get_categorical_variables(new_X)
    categorical <- categorical_list[["categorical"]]
    n_categorical <- categorical_list[["n_categorical"]]
    if (n_categorical > 0) {
      new_X[, categorical] <- sapply(new_X[, categorical], as.numeric)
      new_X[, categorical] <- sapply(new_X[, categorical], function(x) x - 1)
    }
    expanded_z <- cbind(expanded_z, new_X)
    expanded_z <- as.matrix(expanded_z)
  }

  if (fit$criterion == "1se") s <- "lambda.1se"
  if (fit$criterion == "min") s <- "lambda.min"

  if (fit$method == "classif") {
    if (output == "raw") type <- "link"
    if (output == "probability") type <- "response"
    if (output == "class") type <- "class"
  } else if (fit$method == "regr") {
    type <- "link"
  }

  yhat <- glmnet:::predict.cv.glmnet(object = fit$cv_glmnet,
                                    newx = expanded_z, s = s, type = type)

  yhat
}


# betas <- fit[[1]]$alpha[, cvfit$cv[[1]]$i1se]
# test <- second_stage(Z, A, y, betas, method = "classif")
# predict_second_stage(new_Z = Z, fit = test, output = "class")
# test$index
