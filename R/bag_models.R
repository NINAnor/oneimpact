#' Summary of a bag of models
#'
#' Creates a bag of models object from a fitted ensemble of models. This function
#' synthesizes results from multiple model resamples, assigns weights based on validation
#' scores, and computes weighted coefficients and predictions for ensemble modeling.
#'
#' @param fitted `[list]` \cr Object containing fitted models, typically output from
#' `[oneimpact::fit_net_logit()]` or `[oneimpact::fit_net_clogit()]` with multiple resamples.
#' Must contain elements: `$models`, `$metric`, `$formula`, `$method`, `$standardize`, and `$samples`.
#' @param data `[data.frame,tibble]` \cr Complete data set used to fit the models.
#' Used here only to compute summary statistics and covariate scaling information.
#' @param score2weight `[function]` \cr Function to convert validation scores to model weights.
#' Should accept arguments: `x` (model result), `weights_col` (column name), and `score_threshold`.
#' Default is `score2weight_min_invmean`. See `[oneimpact::score2weight_min_invmean()]` and related functions.
#' @param metric `[character]` \cr Name of the evaluation metric to extract from models.
#' Default uses the same metric from the fitted models. Common values: `"AUC"`, `"Cindex"`, `"conditionalAUC"`, 
#' `"conditionalSomersD"`.
#' @param weights_col `[character(1)="validation_score"]` \cr Column name containing validation scores
#' to be converted to weights. One of `"validation_score"` or `"habitat_validation_score"`. The latter
#' is used when habitat-specific validation scores are available, such as in step-selection functions
#' in which one wants to disentangle the effects of movement and habitat.
#' @param score_threshold `[numeric(1)=0.7]` \cr Minimum validation score threshold.
#' Models with validation score below this threshold are assigned zero weight and excluded from
#' the bag. 
#' @param weights_function `[function or NULL]` \cr Function to normalize and transform weights into the
#' range [0,1]. Default is `w_strech_maxmin_squared`. See `[oneimpact::w_strech_maxmin_squared()]` 
#' and related functions.
#'
#' @return A `bag` object (also a `list`) containing:
#' \describe{
#'   \item{n}{Total number of models in the bag}
#'   \item{n_no_errors}{Number of successfully fitted models}
#'   \item{n_errors}{Number of models that failed to fit}
#'   \item{n_above_threshold}{Number of models with weights above zero}
#'   \item{formula}{Model formula used for fitting}
#'   \item{formula_no_strata}{Formula without strata terms}
#'   \item{method}{Modeling method (e.g., "Lasso", "Ridge", "AdaptiveLasso", etc.)}
#'   \item{metric}{Evaluation metric used (e.g. "AUC", "Cindex", etc.)}
#'   \item{weights}{Vector of model weights (normalized)}
#'   \item{coef}{Matrix of raw coefficients for each model}
#'   \item{coef_std}{Matrix of standardized coefficients for each model}
#'   \item{wcoef}{Weighted coefficients}
#'   \item{wcoef_std}{Weighted standardized coefficients}
#'   \item{validation_score}{Matrix of validation scores}
#'   \item{weighted_validation_score}{Weighted validation scores}
#'   \item{covariate_mean_sd}{Data frame with mean and SD of covariates}
#'   \item{data_summary}{Summary statistics of variables in the original data}
#'   \item{errors}{Logical vector indicating which models failed}
#' }
#'
#' @seealso [oneimpact::fit_net_logit()], [oneimpact::fit_net_clogit()],
#' [oneimpact::predict.bag()], [oneimpact::score2weight_min_invmean()],
#' [oneimpact::w_strech_maxmin_squared()]
#'
#' @examples
#' \dontrun{
#' # Assuming you have fitted models from fit_net_logit or fit_net_clogit
#' fitted_models <- bag_fit_net_logit(f = case_ ~ x + y, data = mydata, samples = samples_list)
#' bag <- bag_models(fitted = fitted_models, data = mydata, metric = "AUC")
#' }
#'
#' @name bag_models
#' @export
bag_models <- function(fitted, data,
                       score2weight = score2weight_min_invmean,
                       metric = fitted$metric,
                       weights_col = c("validation_score", "habitat_validation_score")[1],
                       score_threshold = 0.7,
                       weights_function = NULL){

  # function to transform scores to weights, if none is provided
  if (is.null(score2weight)){
    score2weight <- score2weight_min_invmean
  }

  # weighing function for scores
  if (is.null(weights_function)){
    weights_function <- w_strech_maxmin_squared
  }

  # initialize result
  result <- list()

  #lres <- lapply(i, function(i){load(paste0(out_dir, "spat_strat_issf_i", i, ".rda")); return(results)})
  lres <- fitted$models

  # Bag composed of n models
  result$n <- fitted$n

  # formula
  result$formula <- f <- fitted$formula
  wcols <- extract_response_strata(f, covars = TRUE)
  result$formula_no_strata <- as.formula(paste0(wcols$response, " ~ -1+", wcols$covars))

  # method
  # assuming all models follow the same method, must be changed to lapply if not
  result$method <- fitted$method
  # metric
  result$metric <- metric
  # metrics evaluated
  result$metrics_evaluated <- sapply(lres[[1]]$metrics_evaluated, function(x) x$metric)
  # samples
  result$samples <- fitted$samples
  # options
  result$standardize <- fitted$standardize

  # errors
  result$errors <- err <- sapply(lres, function(x) class(x) == "try-error")
  result$error_message <- sapply(1:length(err),
                                 function(i) {
                                   if(err[i]) return(lres[[i]]) else return(NA)
                                 })
  # number of errors
  result$n_errors <- sum(err)
  # final number of valid models
  result$n_no_errors <- result$n - result$n_errors
  if(result$n_no_errors == 1) {
    warning("Only one model in the bag was succefully fitted. The results will be based on this single model.")
  } else {
    if(result$n_no_errors < 1) {
      stop("No models were succefully fitted. Please check and re-fit the models.")
    }
  }

  # parameters for fitting
  result$parms <- lres[[1]]$parms
  # removing redundant elements
  result$parms[c("samples", "i", "metric", "method", "standardize")] <- NULL
  # alpha
  result$alpha <- lres[!err][[1]]$alpha
  # other model outputs
  result$var_names <- lres[!err][[1]]$var_names

  # synthesize results
  if(metric == fitted$metric) {

    # lambda
    lambda <- do.call("cbind", lapply(lres[!err], function(x) { x$lambda } ))
    # coefficients
    coef <- do.call("cbind", lapply(lres[!err], function(x) { x$coef } ))
    coef_std <- do.call("cbind", lapply(lres[!err], function(x) { x$coef_std } ))
    colnames(coef) <- names(lres[!err])
    if(!is.null(coef_std)) colnames(coef_std) <- colnames(coef)

    fit_score <- do.call("cbind", lapply(lres[!err], function(x) { x$train_score} ))
    calibration_score <- do.call("cbind", lapply(lres[!err], function(x) { x$test_score} ))
    validation_score <- do.call("cbind", lapply(lres[!err], function(x) { x$validation_score} ))
    habitat_validation_score <- do.call("cbind", lapply(lres[!err], function(x) { x$habitat_validation_score} ))

  } else {

    # lambda
    lambda <- do.call("cbind", lapply(lres[!err], function(x) { x$lambdas[names(x$lambdas) == metric] } ))
    # coefficients
    coef <- do.call("cbind", lapply(lres[!err], function(x) { x$coefs_all[,colnames(x$coefs_all) == metric, drop = FALSE] } ))
    coef_std <- do.call("cbind", lapply(lres[!err], function(x) { x$coefs_std_all[,colnames(x$coefs_all) == metric, drop = FALSE] } ))
    colnames(coef) <- names(lres[!err])
    if(!is.null(coef_std)) colnames(coef_std) <- colnames(coef)

    fit_score <- do.call("cbind", lapply(lres[!err], function(x) { x$train_score_all[colnames(x$coefs_all) == metric, drop = FALSE] } ))
    calibration_score <- do.call("cbind", lapply(lres[!err], function(x) { x$test_score_all[colnames(x$coefs_all) == metric, drop = FALSE] } ))
    validation_score <- do.call("cbind", lapply(lres[!err], function(x) { x$validation_score_all[,colnames(x$coefs_all) == metric, drop = FALSE] } ))
    habitat_validation_score <- do.call("cbind", lapply(lres[!err], function(x) {
      if(is.null(x$habitat_validation_score_all[[1]])) {
        NULL
      } else {
        x$habitat_validation_score_all[,colnames(x$coefs_all) == metric, drop = FALSE]
      }
    }))

  }

  # weights
  weights <- unlist(lapply(lres[!err], score2weight, col = weights_col, score_threshold = score_threshold))
  weights <- weights_function(weights)
  if(anyNA(weights)) {
    warning("At least one resample has NaN weight. Please check the models' validation scores and the 'score_threshold' parameter.")
  }
  if(all(is.na(weights))) {
    stop("All resamples have NaN weight. Please check the models' validation scores and possibly increase the 'score_threshold' parameter.")
  }

  # lambda
  result$lambda <- lambda
  # weights column
  result$weight_ref <- weights_col
  # weights threshold
  result$weight_threshold <- score_threshold
  # weights
  result$weights <- weights
  # number of models above the threshold
  result$n_above_threshold <- sum(weights > 0)

  # weighted coefficients
  weights_matrix <- matrix(rep(weights, each = nrow(coef)), nrow = nrow(coef))
  wcoef <- coef * as.vector(weights_matrix) # each model
  wcoef_std <- coef_std * as.vector(weights_matrix) # each model

  result$coef <- coef
  result$coef_std <- coef_std
  result$wcoef <- wcoef
  result$wcoef_std <- wcoef_std
  result$fit_score <- fit_score
  result$calibration_score <- calibration_score
  result$validation_score <- validation_score
  result$habitat_validation_score <- habitat_validation_score

  # average validation scores
  if(nrow(result$validation_score) > 1) {
    result$validation_score_summary <- data.frame(min = apply(result$validation_score, 2, min,na.rm = TRUE),
                                                  median = apply(result$validation_score, 2, median,na.rm = TRUE),
                                                  mean = apply(result$validation_score, 2, mean,na.rm = TRUE),
                                                  max = apply(result$validation_score, 2, max,na.rm = TRUE))
  } else {
    result$validation_score_summary <- data.frame(median = apply(result$validation_score,2,mean,na.rm = TRUE))
  }

  #apply(result$validation_score_summary, 2, mean)

  # weighted validation
  result$weighted_validation_score <-  result$validation_score %*% result$weights
  colnames(result$weighted_validation_score) <- "weighted_validation_score"

  if(nrow(result$validation_score) > 1) {

    result$weighted_validation_score_summary <- cbind(
      apply(result$validation_score, 2, min, na.rm = TRUE) %*% result$weights,
      apply(result$validation_score, 2, median,na.rm = TRUE) %*% result$weights,
      apply(result$validation_score, 2, mean,na.rm = TRUE) %*% result$weights,
      apply(result$validation_score, 2, max,na.rm = TRUE) %*% result$weights)
    colnames(result$weighted_validation_score_summary) <- c("min", "median", "mean", "max")
  } else {
    result$weighted_validation_score_summary <-
      apply(result$validation_score,2,min,na.rm = TRUE) %*% result$weights
    colnames(result$weighted_validation_score_summary) <- "mean"
  }

  # covariates summary - based on the whole dataset; different from the ones used for each model in the bag
  all_vars <- all.vars(result$formula_no_strata)
  all_covars <- all_vars[-1]
  # get predictors
  data_covs <- data[, all_covars]
  # select numeric predictors to be standardized
  numeric_covs <- (sapply(data_covs, is.numeric) == TRUE)
  # numeric ones
  data_covs_num <- data_covs[, numeric_covs]
  # standardize
  data_covs_num_std <- lapply(1:ncol(data_covs_num), function(i) scale(data_covs_num[,i]))
  # register mean and sd
  covs_mean_sd <- data.frame(do.call("rbind",lapply(1:length(data_covs_num_std), function(i)
    sapply(c("scaled:center", "scaled:scale"), function(m) attr(data_covs_num_std[[i]], m)))))
  rownames(covs_mean_sd) <- colnames(data_covs_num)
  colnames(covs_mean_sd) <- c("mean", "sd")
  result$covariate_mean_sd <- covs_mean_sd

  # is data numeric?
  classes_numeric <- sapply(data[,all_vars], is.numeric)
  # numeric variables
  data_summary_num <- as.data.frame(apply(na.omit(as.matrix(data[,all_vars[classes_numeric == TRUE]])), 2, data_summary))
  # character variables - use mode
  data_summary_ch <- as.data.frame(apply(na.omit(as.matrix(data[,all_vars[classes_numeric == FALSE]])), 2, data_summary_char))
  names(data_summary_ch) <- all_vars[classes_numeric == FALSE]
  if(nrow(data_summary_ch) > 0) {
    result$data_summary <- cbind(data_summary_num, data_summary_ch)[order(c(which(classes_numeric == TRUE), which(classes_numeric == FALSE)))]
  } else {
    result$data_summary <- data_summary_num
  }

  # identification of numeric covariates
  result$numeric_covs <- lres[[1]]$numeric_covs

  # set class
  class(result) <- c("bag", "list")
  # class(result) <- "bag"

  return(result)
}

#' Compute data summary statistics
#'
#' Utility functions to compute summary statistics for numeric and character variables.
#' Used internally by `bag_models()` to create data summaries.
#'
#' @param x `[numeric or character vector]` \cr Vector of values to summarize.
#' @param na.rm `[logical(1)=TRUE]` \cr Should NA values be removed before computation?
#'
#' @details
#' `data_summary()` computes min, quantiles (1%, 2.5%, 25%, 50%, 75%, 97.5%, 99%),
#' max, mean, and standard deviation for numeric variables.
#'
#' `data_summary_char()` computes mode and quantile-like positions for character/factor
#' variables, returning results with compatible naming.
#'
#' @return Named numeric vector with summary statistics.
#'
#' @keywords internal
#' @name data_summary
#' @export
data_summary <- function(x, na.rm = TRUE){
  y <- c(min(x, na.rm = na.rm), quantile(x, probs=c(0.01, 0.025, 0.25, 0.5, 0.75, 0.975, 0.99), na.rm = na.rm),
         max(x, na.rm = na.rm), mean(x, na.rm = na.rm), sd(x, na.rm = na.rm))
  names(y)[c(1,9,10,11)] <- c("min", "max", "mean", "sd")
  return(y)
}

#' @rdname data_summary
#' @export
data_summary_char <- function(x){
  tab <- table(x)
  tab <- tab[order(tab)]
  nam <- names(tab)
  mode <- nam[which(abs(tab - median(tab)) == min(abs(tab - median(tab))))[1]]
  y <- c(nam[1], rep(mode, 3), mode,
         rep(mode, 3), nam[length(tab)], rep(mode, 2))
  names(y) <- c("min", paste0(c(1, 2.5, 25, 50, 75, 97.5, 99), "%"),
                "max", "mean", "sd")
  return(y)
}

# getmode <- function(v) {
#   uniqv <- unique(v)
#   uniqv[which.max(tabulate(match(v, uniqv)))]
# }

#-----------------
# Weighting functions

#' Helper functions for weight calculation in bag models
#'
#' These functions compute weights for ensemble model combinations based on validation scores.
#' They handle score thresholding and weight normalization differently:
#' - `score2weight_*` functions convert validation scores to weights
#' - `w_strech_*` functions normalize weights after they are computed
#'
#' @param x `[list or numeric]` \cr For `score2weight_*` functions: a model result list with
#' a column specified by `col` parameter. For `w_strech_*` functions: a numeric vector of weights.
#' @param col `[character(1)="validation_score"]` \cr Column name in `x` containing validation scores.
#' @param score_threshold `[numeric(1)=0.7]` \cr Minimum validation score threshold.
#' Models with any validation score below this are assigned zero weight.
#'
#' @details
#' **Score to weight conversion functions:**
#' - `score2weight_mean`: Uses mean of all validation scores
#' - `score2weight_min_mean`: Uses mean only if all scores above threshold
#' - `score2weight_invmean`: Uses inverse of mean (harmonic mean-related)
#' - `score2weight_min_invmean`: Uses inverse mean only if all scores above threshold
#'
#' **Weight stretching/normalization functions:**
#' - `w_strech_maxmin_squared`: Rescales to [0,1], squares, then normalizes to sum to 1
#' - `w_strech_max_squared`: Divides by max, squares, then normalizes to sum to 1
#'
#' @return Numeric vector or value representing normalized weights or summary statistics.
#'
#' @seealso [oneimpact::bag_models()], [oneimpact::data_summary()]
#'
#' @keywords internal
#' @name bag_model_weights
#' @export
w_strech_maxmin_squared <- function(x){
  x <- x-min(x, na.rm = T)
  x <- x/max(x, na.rm = T)
  x <- x^2
  x <- x/sum(x, na.rm = T)
  return(x)
}

#' @rdname bag_model_weights
#' @export
w_strech_max_squared <- function(x){
  # x <- x-min(x, na.rm = T)
  x <- x/max(x, na.rm = T)
  x <- x^2
  x <- x/sum(x, na.rm = T)
  return(x)
}

#' @rdname bag_model_weights
#' @export
score2weight_mean <- function(x, col = "validation_score", score_threshold = 0.7){
  x <- x[[col]]
  x <- mean(x, na.rm = TRUE)
  x <- ifelse(x < score_threshold, 0, x) #set poorly validated models to zero
  return(x)
}

#' @rdname bag_model_weights
#' @export
score2weight_min_mean <- function(x, col = "validation_score", score_threshold = 0.7){
  # x = lres[[14]]
  x <- x[[col]]
  x <- ifelse(any(x < score_threshold), 0, mean(x, na.rm = TRUE)) #set poorly validated models to zero
  return(x)
}

#' @rdname bag_model_weights
#' @export
score2weight_invmean <- function(x, col = "validation_score", score_threshold = 0.7){
  x <- x[[col]]
  x <- 1/mean(1/x, na.rm = TRUE)
  x <- ifelse(x < score_threshold, 0, x) #set poorly validated models to zero
  return(x)
}

#' @rdname bag_model_weights
#' @export
score2weight_min_invmean <- function(x, col = "validation_score", score_threshold = 0.7){
  x <- x[[col]]
  x <- ifelse(any(x < score_threshold), 0, 1/mean(1/x, na.rm = TRUE)) #set poorly validated models to zero
  return(x)
}
