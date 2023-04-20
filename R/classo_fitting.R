#' Perform sparse log-contrast regression
#'
#' Solves the constrained lasso problem using the CLASSO module
#' in Python.  The optimization problem are
#'
#' Regression:
#' minimize_{beta, beta0} 1/(2n) || y - beta0 1_n - Zagg_clr beta ||^2
#'                                     + lamda_max * frac || beta ||_1
#' subject to C beta = 0
#' Classification:
#' minimize_{beta, beta0} max(1 - y_i(beta0 + Z_clr_i * beta), 0)^2 +
#'                               lambda_max * frac || W * beta ||_1
#' subject to C beta = 0
#'
#' Default is C = 1_p^T, but C can be a general matrix.
#'
#' Observe that the tuning parameter is specified through "frac", the fraction
#' of lamda_max (which is the smallest value for which beta is nonzero).
#'
#' @param X n by p matrix containing e.g. log-ratios
#' @param y n vector (response)
#' @param C m by p matrix. Default is a row vector of zeros.
#' @param fraclist (optional) vector of tuning parameter multipliers.
#'    Should be in (0, 1].
#' @param nlam number of tuning parameters (ignored if fraclist non-NULL)
#' @param min_frac smallest value of tuning parameter multiplier (ignored if
#'   fraclist non-NULL)
#' @param w of positive weights of
#'   length ncol(X) (default: all equal to 1).
#' @param method string which estimation method to use should be in
#'   ("regr", "classif", "classif_huber")
#' @param intercept only works for classification! Should the intercept be
#'   fitted. Default is TRUE, set to FALSE if the intercept should not be
#'   included
#' @param normalized if `TRUE` normalize the additional covariates.
#'   In this case the calculation for each covariate / feature:
#'   (X-X_mean) / ||X||_2
#'   The weights will be transformed back to the original scale.
#' @param rho value for huberized classification loss.
#'   Default = 0.0.
#'
#' @export
classo_fitting <- function(X, y, C = NULL, fraclist = NULL,
                            nlam = 20, min_frac = 1e-4,
                            method = c("regr", "classif", "classif_huber"),
                            w = NULL, intercept = TRUE,
                            normalized = TRUE, rho = 0.0) {
  n <- length(y)
  stopifnot(nrow(X) == n)
  p <- ncol(X)
  if (is.null(C)) C <- matrix(0, nrow = 1, ncol = p)

  # check method
  method <- match.arg(method)
  method_check <- check_method(method = method, y = y, rho = rho)

  classification <- method_check$classification
  y <- method_check$y


  # normalize data
    if (normalized) {
      # call the normalization helper function
      normalized_values <-
        normalization_additional_covariates(additional_covariates =
                                            as.data.frame(X),
                                            p_x = p,
                                            intercept = intercept)
      X <- normalized_values$X
    } else {
      # get the number of categorical variables if no normalization is applied
      categorical_list <- get_categorical_variables(as.data.frame(X))
      categorical <- categorical_list[["categorical"]]
      n_categorical <- categorical_list[["n_categorical"]]
      if (n_categorical > 0) {
        additional_covariates[, categorical] <-
          transform_categorical_variables(as.data.frame(X), categorical)
    }
    # define the weights as 1 for every non compositional covariates
    # since c-lasso can now handle weights and we do not need to transform
    # them anymore and pass them directly to the solver
    # the weights for the compositional effects is taken into account by
    # modifying C and the weights for the non-compositional effects through
    # w in c-lasso
    }


  # CLASSO can solve a problem of the form
  #
  # (based on https://github.com/Leo-Simpson/c-lasso)
  # [R1] Standard constrained Lasso regression:
  # minimize_{beta} || y - X beta ||^2 + lam ||beta||_1
  # subject to C beta = 0
  #
  # Define
  # yt = y - ybar 1_n
  # M = Z_clr - 1_n v^T, where v = colMeans(Z_clr)
  #
  if (is.null(fraclist))
    fraclist <- exp(seq(0, log(min_frac), length = nlam))


  ybar <- mean(y)
  if (classification) {
    # for classification we do not need to scale the outcome
    yt <- y
  } else {
    # scale y
    # ybar <- mean(y)
    # yt <- y - ybar
  }

  if (!classification) intercept <- TRUE

  fit <- list()
  X_classo <- as.matrix(X)

  # set up CLASSO problem:
  prob <- classo$classo_problem(X = X_classo,
                                C = C,
                                y = array(yt))
  prob$formulation$classification <- classification
  prob$formulation$concomitant <- FALSE
  if (classification & intercept) {
    prob$formulation$intercept <- TRUE
  }
  if (method == "classif_huber") {
    prob$formulation$huber <- TRUE
    prob$formulation$rho_classification <- rho
  } else {
    prob$formulation$huber <- FALSE
  }
  prob$model_selection$PATH <- TRUE
  prob$model_selection$CV <- FALSE
  prob$model_selection$LAMfixed <- FALSE
  prob$model_selection$StabSel <- FALSE
  prob$model_selection$PATHparameters$lambdas <- fraclist
  if (!is.null(w)) prob$formulation$w <- w


  # solve  it
  prob$solve()
  # extract outputs

  beta <- (prob$solution$PATH$BETAS)
  # purrr::map(prob$solution$PATH$BETAS, as.numeric) %>%
  #   rlang::set_names(paste0("V", 1:length(prob$solution$PATH$BETAS))) %>%
  #   dplyr::bind_cols() %>%
  #   as.matrix()
  if (intercept) {
    # c-lasso can estimate beta0 --> select first column (estimated beta0)
    # delete the first column afterwards
    beta0 <- beta[, 1]
    beta <- beta[, -1]
  }
  beta <- t(beta)
  lambda_classo <- prob$model_selection$PATHparameters$lambdas
#  if (!classification) beta0 <- ybar - crossprod(beta[1:p, ], v)
  if (!intercept) beta0 <- rep(0, times = length(lambda_classo))



  if (normalized && (normalized_values$n_numeric != 0)) {
      # rescale betas for numerical values
      # rescale only if beta not 0
      beta <- rescale_betas(
        beta = beta,
        p_x = p,
        p = 0,
        n_numeric = normalized_values$n_numeric,
        categorical = normalized_values$categorical,
        xs = normalized_values$xs,
        xm = normalized_values$xm
      )
      rownames(beta)[1:p] <- colnames(X)
  } else {
    rownames(beta)[1:p] <- colnames(X)
  }

  list(beta0 = beta0,
       beta = beta,
       C = C,
       fraclist = lambda_classo, # / (2 * n),
       fit_classo = prob,
       refit = FALSE,
       w = w,
       method = method,
       intercept = intercept,
       rho = rho,
       normalized = normalized)
}
