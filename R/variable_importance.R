#' Computes variable importance
#'
#' Computation of variable importance from a bag of models. Variable importance can
#' be computed either by dropping variables or through permutation of the variables values
#' (parameter `type`), typically by evaluating the effects on the model evaluation metric
#' in the validation set(s). If `type = "drop"`, each variable is dropped from the model
#' at a time and the variation in model evaluation metric is computed. If `type = "permutation"`,
#' The observations of each variable are permutated and the variation in model evaluation
#' metric is computed.
#'
#' @param x `[list]` \cr Bag of models, result of [oneimpact::bag_models()]. It contains
#' multiple information from the models, such as formula, weights, coefficients,
#' and metric used to evaluate the models.
#' @param data `[data.frame]` \cr Complete data set to which the models were applied.
#' @param type `[character(1)="drop"]{"drop", "permutation"}` \cr Type of computation
#' for variable importance. If `type = "drop"` (default), each variable is dropped from
#' the model at a time and the variation in model evaluation metric is computed.
#' If `type = "permutation"`, the observations of each variable are permutated and the
#' variation in the model evaluation metric is computed.
#' @param samples `[list]` \cr List of samples used to fit the models in the bag.
#' The list contains at least three elements: train, test,
#' and validate. Each elements might have several elements, each representing
#' the lines of `data` to be sampled for each resample. Typically, this is computed by
#' the function [oneimpact::create_resamples()].
#' @param n_permutations `[numeric(1)=100]` \cr Number of permutations, if
#' `type = "permutation"`.
#' @param order `[character,logical(1)="desc"]{"desc", "asc", FALSE}` \cr Whether or
#' not to order the output variables according to descending (`order = "desc"`) or
#' ascending order of variable importance (`order = "asc"`). If `FALSE`, the variables
#' are shown in the same order as present in the bag of models, `x`.
#' @param plot `[logical(1)=FALSE]` \cr Should variable importance be plotted? Default
#' is `FALSE`.
#' @param remove_threshold `[numeric(1)]` \cr Threshold for excluding variable with
#' little importance in the variable importance plot (i.e. only considered if `plot = TRUE`).
#' See more in [oneimpact::plot_importance()].
#'
#' @seealso For plotting variable importance, see [oneimpact::plot_importance()].
#' @export
variable_importance <- function(x,
                                data,
                                samples = NULL,
                                type = c("drop", "permutation")[1],
                                n_permutations = 100,
                                order = c("desc", "asc", FALSE)[1],
                                variable_block = NULL,
                                metric = NULL,
                                plot = FALSE,
                                remove_threshold = 0) {
  # get info
  f <- x$formula
  wghts <- x$weights
  coefs <- x$coef
  if(is.null(metric)) metric <- x$metric
  # set sample(s)
  ss <- 1

  # case
  case <- extract_response_strata(f)$response
  # strata?
  strat <- extract_response_strata(f)$strata

  # should training data be used?
  if(!is.null(samples)) {
    if(strat == "") {
      data <- data[samples$validate[[ss]],]
    } else {
      data <- data[data[[strat]] %in% data[data[[case]] == 1,][[strat]][samples$validate[[ss]]],]
    }
    ## here we select the first sample; it could be for all, and an average over samples
    # or for all data altogether, with not samples
    ### change here for conditional logistic regression
  }

  if(!is.null(variable_block)) {
    # check length of variable block
  }

  # model matrix, prediction, y and strata
  mm <- model.matrix(x$formula_no_strata, data)
  pred <- as.vector((mm %*% coefs) %*% wghts)
  y <- data[, case]

  if(strat == "") {
    strat <- 1
  } else {
    strat <- data[,strat]
  }

  baseline <- metric(data.frame(x = pred, y, strat))

  if (type == "permutation"){
    if(is.null(variable_block)) {
      test <- unlist(lapply(c(1:ncol(mm)), permutated_concordance, mm = mm, y = y, strat = strat,
                            coefs = coefs, wghts = wghts, n_permutations = n_permutations,
                            baseline_pred = pred, metric = metric))
    } else {
      test <- unlist(lapply(c(1:ncol(mm)), permutated_concordance, mm = mm, y = y, strat = strat,
                            coefs = coefs, wghts = wghts, n_permutations = n_permutations,
                            baseline_pred = pred, metric = metric, variable_block = variable_block))
    }
  } else {
    if(is.null(variable_block)) {
      test <- unlist(lapply(c(1:ncol(mm)), dropped_concordance, mm = mm, y = y, strat = strat,
                            coefs = coefs, wghts = wghts, baseline_pred = pred, metric = metric))
    } else {
      test <- unlist(lapply(c(1:length(unique(variable_block))),
                            dropped_concordance, mm = mm, y = y, strat = strat,
                            coefs = coefs, wghts = wghts, baseline_pred = pred, metric = metric,
                            variable_block = variable_block))
    }

  }
  test <- ifelse(test > baseline, baseline, test)
  test <- baseline - test
  test <- test/baseline
  test <- test/sum(test)

  if(is.null(variable_block)) {
    names(test) <- colnames(mm)
  } else {
    names(test) <- unique(variable_block)
  }

  # order?
  if(order == "desc") {
    test <- test[order(-test)]
  } else {
    if(order == "asc") {
      test <- test[order(test)]
    }
  }

  # plot?
  if(plot) print(plot_importance(test, remove_threshold = remove_threshold))

  return(test)
}

# add variable block here??
permutated_concordance <- function(i, mm, y, strat, coefs, wghts,
                                   n_permutations, baseline_pred,
                                   metric, variable_block = NULL){
  tmp <- unlist(lapply(c(1:n_permutations),
                       function(x, i, mm, y, strat, coefs, wghts, baseline_pred){
                         # matrix n_obs rows x r_resamples cols
                         mat_coef <- (matrix(rep(coefs[i,], each = nrow(mm)), nrow = nrow(mm)))
                         # change in prediction value
                         pred <- baseline_pred + (((mm[sample.int(nrow(mm)),i]-mm[,i]) * mat_coef) %*% wghts)
                         return(metric(data.frame(x=pred, y, strat), warnings=F))},
                       i=i, mm=mm, y=y, strat=strat, coefs=coefs, wghts=wghts, baseline_pred=baseline_pred))
  return(mean(tmp))
}

permutated_concordance_old <- function(i, mm, y, strat, coefs, wghts,
                                       n_permutations, baseline_pred,
                                       metric, variable_block = NULL){
  tmp <- unlist(lapply(c(1:n_permutations),
                       function(x, i, mm, y, strat, coefs, wghts, baseline_pred){
                         mm[,i] <- mm[sample.int(nrow(mm)),i];
                         pred <- as.vector((mm %*% coefs) %*% wghts);
                         return(metric(data.frame(x=pred, y, strat), warnings=F))},
                       i=i, mm=mm, y=y, strat=strat, coefs=coefs, wghts=wghts, baseline_pred=baseline_pred))
  return(mean(tmp))
}
#
# permutated_concordance(7, mm, y, strat, coefs, wghts, n_permutations=1, baseline_pred, metric)
# permutated_concordance_old(7, mm, y, strat, coefs, wghts, n_permutations=1, metric)#
# microbenchmark::microbenchmark(permutated_concordance(7, mm, y, strat, coefs, wghts, n_permutations=1, baseline_pred, metric),
#                                permutated_concordance_old(7, mm, y, strat, coefs, wghts, n_permutations=1, baseline_pred, metric))

dropped_concordance_new <- function(i, mm, y, strat, coefs, wghts,
                                    n_permutations, baseline_pred,
                                    metric, variable_block = NULL){

  # matrix n_obs rows x r_resamples cols
  weights_nonzero <- which(wghts > 0)

  if(!is.null(variable_block)) {
    ### NEED TO CHECK THIS HERE
    mat_coef <- (matrix(rep(coefs[variable_block == unique(variable_block)[i], weights_nonzero], each = nrow(mm)), nrow = nrow(mm)))
    pred <- baseline_pred - ((mm[, variable_block == unique(variable_block)[i]] * mat_coef) %*% wghts[weights_nonzero])
  } else {
    mat_coef <- (matrix(rep(coefs[i,weights_nonzero], each = nrow(mm)), nrow = nrow(mm)))
    pred <- baseline_pred - ((mm[, i] * mat_coef) %*% wghts[weights_nonzero])
  }

  return(metric(data.frame(x=pred, y, strat), warnings=F))
}

dropped_concordance <- function(i, mm, y, strat, coefs, wghts,
                                n_permutations, baseline_pred,
                                metric, variable_block = NULL){
  if(!is.null(variable_block)) {
    mm[,variable_block == unique(variable_block)[i]] <- 0
  } else {
    mm[,i] <- 0
  }

  weights_nonzero <- which(wghts > 0)

  pred <- as.vector((mm %*% coefs[,weights_nonzero]) %*% wghts[weights_nonzero]);
  return(metric(data.frame(x=pred, y, strat), warnings=F))
}

# dropped_concordance(7, mm, y, strat, coefs, wghts, n_permutations=1, baseline_pred, metric)
# dropped_concordance_old(7, mm, y, strat, coefs, wghts, n_permutations=1, metric)#
# microbenchmark::microbenchmark(dropped_concordance(7, mm, y, strat, coefs, wghts, n_permutations=1, baseline_pred, metric),
#                                dropped_concordance_old(7, mm, y, strat, coefs, wghts, n_permutations=1, metric),
#                                check = "identical")


# # not used now, replaced by dropped_concordance which allows for other metrics
# permutatedBoyce <- function(i, mm, y, strat, coefs, wghts, n_permutations){
#   tmp <- unlist(lapply(c(1:n_permutations),
#                        function(x, i, mm, y, strat, coefs, wghts){
#                          mm[,i] <- mm[sample.int(nrow(mm)),i];
#                          pred <- as.vector((mm %*% coefs) %*% wghts);
#                          return(conditionalBoyce(data.frame(x=pred, y, strat), warnings=F))}, i=i, mm=mm, y=y, strat=strat, coefs=coefs, wghts=wghts))
#   return(mean(tmp))
# }
# # not used now, replaced by dropped_concordance which allows for other metrics
# droppedBoyce <- function(i, mm, y, strat, coefs, wghts, n_permutations){
#   mm[,i] <- 0
#   pred <- as.vector((mm %*% coefs) %*% wghts);
#   return(conditionalBoyce(data.frame(x=pred, y, strat), warnings=F))
# }


#' Plot variable importance
#'
#' @param importance `[numeric]` \cr Numeric vector showing the importance score of
#' each variable, resulting from [oneimpact::variable_importance()].
#' @param remove_threshold `[numeric(1)=0]` \cr Only importance scores above this
#' level will be plotted. Default is 0, in case which all null variables will be
#' removed. To plot all variables, set remove_threshold as -1 or any other negative value.
#' @param normalize `[logical(1)=TRUE]` \cr If `TRUE`, the variable importance scores are
#' dividing the the maximum score, so that the score of the most important variable
#' is set to 1.
#'
#' @export
plot_importance <- function(importance, remove_threshold = 0, normalize = TRUE) {

  # normalization
  if(normalize)
    imp <- importance/max(importance) else
      imp <- importance

  # data frame
  df <- data.frame(var = factor(names(importance), levels = names(importance), ordered = TRUE),
                   importance = imp)
  df <- df[df$importance > remove_threshold,]

  # plot
  p <- ggplot2::ggplot(data = df, ggplot2::aes(x = var, y = importance)) +
    ggplot2::geom_bar(stat="identity", fill="steelblue") +
    # scale_y_continuous(trans='log10') +
    ggplot2::theme_minimal() +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Variable", y = "Importance")
  p

}
