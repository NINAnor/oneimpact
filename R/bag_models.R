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
                       score2weight = score2weight_invmean,
                       weights_col = c("validation_score", "habitat_validation_score")[1],
                       score_threshold = 0.7,
                       weights_function = NULL,
                       out_dir){

  # function to transform scores to weights, if none is provided
  if (is.null(score2weight)){
    score2weight <- score2weight_invmean
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
  wcols <- extract_response_strata(f, other_vars = TRUE)
  result$mm_formula <- as.formula(paste0(wcols$response, " ~ -1+", wcols$other_vars))

  # method
  # assuming all models follow the same method, must be changed to lapply if not
  result$method <- fitted$method
  # weights
  result$weight_ref <- weights_col
  # weight threshold
  result$weight_threshold <- score_threshold
  # metric
  result$metric <- fitted$metric

  # errors
  result$errors <- err <- sapply(lres, function(x) class(x) == "try-error")
  result$error_message <- sapply(1:length(err),
                                 function(i) {
                                   if(err[i]) return(lres[[i]]) else return(NA)
                                 })
  # number of errors
  result$n_errors <- sum(err)
  # final number of valid models
  result$n_final <- result$n - result$n_errors
  if(result$n_final == 1) {
    warning("Only one model in the bag was succefully fitted. The results will be based on this single model.")
  } else {
    if(result$n_final < 1) {
      stop("No models were succefully fitted. Please check and re-fit the models.")
    }
  }

  # should we unstandardize the coefs here!?!

  # synthesize results
  result$coef <- do.call("cbind", lapply(lres[!err], function(x) { x$coef } ))
  colnames(result$coef) <- names(lres[!err])
  result$var_names <- lres[!err][[1]]$var_names
  result$fit_score <- do.call("cbind", lapply(lres[!err], function(x) { x$fit_score} ))
  result$calibration_score <- do.call("cbind", lapply(lres[!err], function(x) { x$calibration_score} ))
  result$validation_score <- do.call("cbind", lapply(lres[!err], function(x) { x$validation_score} ))
  result$habitat_validation_score <- do.call("cbind", lapply(lres[!err], function(x) { x$habitat_validation_score} ))

  # weights
  result$weights <- unlist(lapply(lres[!err], score2weight, col = weights_col, score_threshold = score_threshold))
  result$weights <- weights_function(result$weights)
  if(anyNA(result$weights)) {
    warning("At least one resample has NaN weight. Please check the models' validation scores the the 'score_threshold' parameter.")
  }

  # average validation scores
  if(nrow(result$validation_score) > 1) {
    result$validation_score_summary <- data.frame(min = apply(result$validation_score,2,min,na.rm = TRUE),
                                                  median = apply(result$validation_score,2,median,na.rm = TRUE),
                                                  mean = apply(result$validation_score,2,mean,na.rm = TRUE),
                                                  max = apply(result$validation_score,2,max,na.rm = TRUE))
  } else {
    result$validation_score_summary <- data.frame(median = apply(result$validation_score,2,mean,na.rm = TRUE))
  }

  #apply(result$validation_score_summary, 2, mean)

  # weighted validation
  if(nrow(result$validation_score) > 1) {
    result$avg_weighted_validation_score <- cbind(
      apply(result$validation_score,2,min, na.rm = TRUE) %*% result$weights,
      apply(result$validation_score,2,median,na.rm = TRUE) %*% result$weights,
      apply(result$validation_score,2,mean,na.rm = TRUE) %*% result$weights,
      apply(result$validation_score,2,max,na.rm = TRUE) %*% result$weights)
    colnames(result$avg_weighted_validation_score) <- c("min", "median", "mean", "max")
  } else {
    result$avg_weighted_validation_score <-
      apply(result$validation_score,2,min,na.rm = TRUE) %*% result$weights
    colnames(result$avg_weighted_validation_score) <- "mean"
  }

  # covariates summary
  all_vars <- all.vars(result$mm_formula)
  classes <- sapply(data[,all_vars], class)
  # numeric variables
  data_summary_num <- as.data.frame(apply(na.omit(as.matrix(data[,all.vars(result$mm_formula)[classes == "numeric"]])), 2, data_summary))
  # character variables - use mode
  data_summary_ch <- as.data.frame(apply(na.omit(as.matrix(data[,all.vars(result$mm_formula)[classes != "numeric"]])), 2, data_summary_char))
  names(data_summary_ch) <- all_vars[classes != "numeric"]
  result$data_summary <- cbind(data_summary_num, data_summary_ch)[order(c(which(classes == "numeric"), which(classes != "numeric")))]

  # identification of numeric covariates
  result$numeric_covs <- fitted$numeric_covs

  # set class
  class(result) <- c("bag", "list")
  # class(result) <- "bag"

  return(result)
}

data_summary <- function(x){
  y <- c(min(x), quantile(x, probs=c(0.01, 0.025, 0.25, 0.5, 0.75, 0.975, 0.99)), max(x), mean(x), sd(x))
  names(y)[c(1,9,10,11)] <- c("min", "max", "mean", "sd")
  return(y)
}

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
  x <- x[[weights_col]]
  x <- mean(x, na.rm = TRUE)
  x <- ifelse(x < score_threshold, 0, x) #set poorly validated models to zero
  return(x)
}

#' @rdname bag_models
#' @export
score2weight_invmean <- function(x, col = "validation_score", score_threshold = 0.7){
  x <- x[[weights_col]]
  x <- 1/mean(1/x, na.rm = TRUE)
  x <- ifelse(x < score_threshold, 0, x) #set poorly validated models to zero
  return(x)
}
