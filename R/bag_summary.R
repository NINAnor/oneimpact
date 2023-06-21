#' Summary of a bag of models
#'
#' @param f Formula (full).
#' @param data Input data (whole dataset)
#' @param score2weight Function to set validation scores into weights, with two arguments:
#' x, the result of one model of the bag, and col, the column to be used for setting the
#' scores. See the argument `weights_col`.
#' @param weights_col Column to use for scores. One of `"validation_score"` or
#' `"habitat_validation_score"`.
#'
#' @export
bag_summary <- function(i, f, data,
                        score2weight = NULL,
                        weights_col = c("validation_score", "habitat_validation_score")[1],
                        score_threshold = 0.7,
                        weights_function = NULL, out_dir){

  # function to transform scores to weights, if none is provided
  if (is.null(score2weight)){
    score2weight <- function(x, col = weights_col){
      x <- x[[weights_col]]
      x <- 1/mean(1/x)
      x <- ifelse(score_threshold < 0.7, 0, x) #set poorly validated models to zero
      return(x)
    }
  }

  # weighing function for scores
  if (is.null(weights_function)){
    weights_function <- function(x){
      x <- x-min(x)
      x <- x/max(x)
      x <- x^2
      x <- x/sum(x)
      return(x)
    }
  }

  result <- list()

  #lres <- lapply(i, function(i){load(paste0(out_dir, "spat_strat_issf_i", i, ".rda")); return(results)})
  lres <- fittedl
  result$coef <- do.call("cbind", lapply(lres, function(x) { x$coef } ))
  result$var_names <- lres[[1]]$var_names
  result$validation_score <- do.call("cbind", lapply(lres, function(x) { x$validation_score} ))
  result$habitat_validation_score <- do.call("cbind", lapply(lres, function(x) { x$habitat_validation_score} ))
  result$weights <- unlist(lapply(lres, score2weight, col = weights_col))
  result$weights <- weights_function(result$weights)

  result$formula <- f
  wcols <- extract_response_strata(f, other_vars = TRUE)
  result$mm_formula <- as.formula(paste0(wcols$case, " ~ -1+", wcols$other_vars))

  result$data_summary <- as.data.frame(apply(na.omit(data[,all.vars(result$mm_formula)]), 2, data_summary))

  return(result)
}

data_summary <- function(x){
  y <- c(min(x), quantile(x, probs=c(0.01, 0.025, 0.25, 0.5, 0.75, 0.975, 0.99)), max(x), mean(x), sd(x))
  names(y)[c(1,9,10,11)] <- c("min", "max", "mean", "sd")
  return(y)
}
