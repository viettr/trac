#' Perform cross validation for tuning parameter selection for sparse log contrast
#'
#' This function is to be called after calling \code{\link{sparse_log_contrast}}.
#' It performs \code{nfold}-fold cross validation.
#'
#' @param fit output of \code{\link{sparse_log_contrast}} function.
#' @param X,y same arguments as passed to
#'   \code{\link{sparse_log_contrast}}. C will be taken from fit object.
#' @param folds a partition of \code{1:nrow(Z)}.
#' @param nfolds number of folds for cross-validation
#' @param summary_function how to combine the errors calculated on each
#' observation within a fold (e.g. mean or median)
#' @export
cv_classo_fitting <- function(fit, X, y, folds = NULL, nfolds = 5,
                               summary_function = stats::median) {
  n <- nrow(X)
  p <- ncol(X)
  stopifnot(length(y) == n)
  if (is.null(folds)) folds <- ggb:::make_folds(n, nfolds)
  else
    nfolds <- length(folds)
  cv <- list()
  fit_folds <- list() # save this to reuse by log-ratio's cv function
  errs <- matrix(NA, ncol(fit$beta), nfolds)
  predicted_values <- matrix(NA, 0, ncol(fit$beta))

  for (i in seq(nfolds)) {
    cat("fold", i, fill = TRUE)
    # add for backward compatibility
    if (is.null(fit$method)) fit$method <- "regr"
    if (is.null(fit$w)) {
      fit$w <- NULL
    }
    if (is.null(fit$rho)) fit$rho <- 0
    if (is.null(fit$normalized)) fit$normalized <- FALSE
    # train on all but i-th fold (and use settings from fit):
    fit_folds[[i]] <- classo_fitting(X[-folds[[i]], ],
                                      y[-folds[[i]]],
                                      fit$C,
                                      fraclist = fit$fraclist,
                                      w = fit$w,
                                      method = fit$method,
                                      rho = fit$rho,
                                      normalized = fit$normalized)
    if (fit$refit) stop("Not yet supported.")
    if (fit$method == "regr" | is.null(fit$method)) {
      errs[, i] <- apply((predict_trac(
        list(fit_folds[[i]]),
        X[folds[[i]], ])[[1]] - y[folds[[i]]])^2,
        2, summary_function
      )
    }

    if (fit$method == "classif" |
        fit$method == "classif_huber") {
      # loss: max(0, 1 - y_hat * y)^2
      er <- sign(predict_trac(list(fit_folds[[i]]),
                              X[folds[[i]],])[[1]]) !=
        c(y[folds[[i]]])
      errs[, i] <- colMeans(er)
    }

    predicted_values <- rbind(predicted_values,
                              predict_trac(
                                list(fit_folds[[i]]),
                                X[folds[[i]], ])[[1]])
  }
  m <- rowMeans(errs)
  se <- apply(errs, 1, stats::sd) / sqrt(nfolds)
  ibest <- which.min(m)
  i1se <- min(which(m < m[ibest] + se[ibest]))
  cv <- list(errs = errs, m = m, se = se,
             lambda_best = fit$fraclist[ibest], ibest = ibest,
             lambda_1se = fit$fraclist[i1se], i1se = i1se,
             fraclist = fit$fraclist,
             nonzeros = colSums(abs(fit$beta) > 1e-5),
             fit_folds = fit_folds,
             predicted_values = predicted_values
  )
  list(cv = cv, folds = folds)
}
