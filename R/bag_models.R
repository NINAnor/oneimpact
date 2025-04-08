#' Summary of a bag of models
#'
#' ### Add option to read the models from a set of rds/rda files.
#'
#' @param data Input data (whole dataset)
#' @param score2weight Function to set validation scores into weights, with two arguments:
#' x, the result of one model of the bag, and col, the column to be used for setting the
#' scores. See the argument `weights_col`.
#' @param weights_col Column to use for scores. One of `"validation_score"` or
#' `"habitat_validation_score"`.
#'
#' @name bag_models
#' @export
bag_models <- function(fitted, data,
                       score2weight = score2weight_min_invmean,
                       metric = fitted$metric,
                       weights_col = c("validation_score", "habitat_validation_score")[1],
                       score_threshold = 0.7,
                       weights_function = NULL,
                       out_dir){

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
  # metric
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
  # other model outputs
  result$var_names <- lres[!err][[1]]$var_names

  # synthesize results
  if(metric == fitted$metric) {

    # coefficients
    coef <- do.call("cbind", lapply(lres[!err], function(x) { x$coef } ))
    coef_std <- do.call("cbind", lapply(lres[!err], function(x) { x$coef_std } ))
    colnames(coef) <- names(lres[!err])
    if(!is.null(coef_std)) colnames(coef_std) <- colnames(coef)

    fit_score <- do.call("cbind", lapply(lres[!err], function(x) { x$fit_score} ))
    calibration_score <- do.call("cbind", lapply(lres[!err], function(x) { x$calibration_score} ))
    validation_score <- do.call("cbind", lapply(lres[!err], function(x) { x$validation_score} ))
    habitat_validation_score <- do.call("cbind", lapply(lres[!err], function(x) { x$habitat_validation_score} ))

  } else {

    # coefficients
    coef <- do.call("cbind", lapply(lres[!err], function(x) { x$coefs_all[,colnames(x$coefs_all) == metric, drop = FALSE] } ))
    coef_std <- do.call("cbind", lapply(lres[!err], function(x) { x$coefs_std_all[,colnames(x$coefs_all) == metric, drop = FALSE] } ))
    colnames(coef) <- names(lres[!err])
    if(!is.null(coef_std)) colnames(coef_std) <- colnames(coef)

    fit_score <- do.call("cbind", lapply(lres[!err], function(x) { x$train_score_all[colnames(x$coefs_all) == metric, drop = FALSE] } ))
    calibration_score <- do.call("cbind", lapply(lres[!err], function(x) { x$test_score_all[colnames(x$coefs_all) == metric, drop = FALSE] } ))
    validation_score <- do.call("cbind", lapply(lres[!err], function(x) { x$validation_score_all[,colnames(x$coefs_all) == metric, drop = FALSE] } ))
    habitat_validation_score <- do.call("cbind", lapply(lres[!err], function(x) { x$habitat_validation_score_all[,colnames(x$coefs_all) == metric, drop = FALSE] } ))

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

#' @export
data_summary <- function(x, na.rm = TRUE){
  y <- c(min(x, na.rm = na.rm), quantile(x, probs=c(0.01, 0.025, 0.25, 0.5, 0.75, 0.975, 0.99), na.rm = na.rm),
         max(x, na.rm = na.rm), mean(x, na.rm = na.rm), sd(x, na.rm = na.rm))
  names(y)[c(1,9,10,11)] <- c("min", "max", "mean", "sd")
  return(y)
}

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

#' @rdname bag_models
#' @export
w_strech_maxmin_squared <- function(x){
  x <- x-min(x, na.rm = T)
  x <- x/max(x, na.rm = T)
  x <- x^2
  x <- x/sum(x, na.rm = T)
  return(x)
}

#' @rdname bag_models
#' @export
w_strech_max_squared <- function(x){
  # x <- x-min(x, na.rm = T)
  x <- x/max(x, na.rm = T)
  x <- x^2
  x <- x/sum(x, na.rm = T)
  return(x)
}

#' @rdname bag_models
#' @export
score2weight_mean <- function(x, col = "validation_score", score_threshold = 0.7){
  x <- x[[col]]
  x <- mean(x, na.rm = TRUE)
  x <- ifelse(x < score_threshold, 0, x) #set poorly validated models to zero
  return(x)
}

#' @rdname bag_models
#' @export
score2weight_min_mean <- function(x, col = "validation_score", score_threshold = 0.7){
  # x = lres[[14]]
  x <- x[[col]]
  x <- ifelse(any(x < score_threshold), 0, mean(x, na.rm = TRUE)) #set poorly validated models to zero
  return(x)
}

#' @rdname bag_models
#' @export
score2weight_invmean <- function(x, col = "validation_score", score_threshold = 0.7){
  x <- x[[col]]
  x <- 1/mean(1/x, na.rm = TRUE)
  x <- ifelse(x < score_threshold, 0, x) #set poorly validated models to zero
  return(x)
}

#' @rdname bag_models
#' @export
score2weight_min_invmean <- function(x, col = "validation_score", score_threshold = 0.7){
  x <- x[[col]]
  x <- ifelse(any(x < score_threshold), 0, 1/mean(1/x, na.rm = TRUE)) #set poorly validated models to zero
  return(x)
}
