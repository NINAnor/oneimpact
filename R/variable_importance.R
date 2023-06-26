#' Computes variable importance
#'
#' @param x `list` \cr Bag of models, result of [oneimpact::bag_models()]. It contains
#' multiple information from the models, such as formula, weights, coefficients,
#' and metric used to evaluate the models.
#'
#' @export
variable_importance <- function(x,
                                data,
                                type = c("drop", "permutation")[1],
                                nb_permutations = 100,
                                order = c("desc", "asc", FALSE)[1],
                                plot = TRUE,
                                remove_threshold = 0) {
  # get info
  f <- x$formula
  wghts <- x$weights
  coefs <- x$coef
  metric <- x$metric

  # model matrix, prediction, y and strata
  mm <- model.matrix(x$mm_formula, data)
  pred <- as.vector((mm %*% coefs) %*% wghts)
  y <- data[,extract_response_strata(f)$response]

  strat <- extract_response_strata(f)$strata
  if(strat == "") {
    strat <- 1
  } else {
    strat <- data[,strat]
  }

  baseline <- metric(data.frame(x=pred, y, strat))

  if (type=="permutation"){
    test <- unlist(lapply(c(1:ncol(mm)), permutated_concordance, mm=mm, y=y, strat=strat,
                          coefs=coefs, wghts=wghts, nb_permutations=nb_permutations,
                          metric = metric))
  }else{
    test <- unlist(lapply(c(1:ncol(mm)), dropped_concordance, mm=mm, y=y, strat=strat,
                          coefs=coefs, wghts=wghts, metric = metric))
  }
  test <- ifelse(test>baseline, baseline, test)
  test <- baseline - test
  test <- test/baseline
  test <- test/sum(test)
  names(test) <- colnames(mm)

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

permutatedBoyce <- function(i, mm, y, strat, coefs, wghts, nb_permutations){
  tmp <- unlist(lapply(c(1:nb_permutations),
                       function(x, i, mm, y, strat, coefs, wghts){
                         mm[,i] <- mm[sample.int(nrow(mm)),i];
                         pred <- as.vector((mm %*% coefs) %*% wghts);
                         return(conditionalBoyce(data.frame(x=pred, y, strat), warnings=F))}, i=i, mm=mm, y=y, strat=strat, coefs=coefs, wghts=wghts))
  return(mean(tmp))
}

permutated_concordance <- function(i, mm, y, strat, coefs, wghts, nb_permutations, metric){
  tmp <- unlist(lapply(c(1:nb_permutations),
                       function(x, i, mm, y, strat, coefs, wghts){
                         mm[,i] <- mm[sample.int(nrow(mm)),i];
                         pred <- as.vector((mm %*% coefs) %*% wghts);
                         return(metric(data.frame(x=pred, y, strat), warnings=F))}, i=i, mm=mm, y=y, strat=strat, coefs=coefs, wghts=wghts))
  return(mean(tmp))
}

droppedBoyce <- function(i, mm, y, strat, coefs, wghts, nb_permutations){
  mm[,i] <- 0
  pred <- as.vector((mm %*% coefs) %*% wghts);
  return(conditionalBoyce(data.frame(x=pred, y, strat), warnings=F))
}

dropped_concordance <- function(i, mm, y, strat, coefs, wghts, nb_permutations, metric){
  mm[,i] <- 0
  pred <- as.vector((mm %*% coefs) %*% wghts);
  return(metric(data.frame(x=pred, y, strat), warnings=F))
}

#' Plot variable importance
#'
#' @param importance `[numeric]` \cr Numeric vector showing the importance score of
#' each variable, resulting from [oneimpact::variable_importance()].
#' @param remove_threshold `[numeric(1)=0]` \cr Only importance scores above this
#' level will be plotted. Default is 0, in case which all null variables will be
#' removed. To plot all variables, set remove_threshold as -1 or any other negative value.
#'
#' @export
plot_importance <- function(importance, remove_threshold = 0) {

  # data frame
  df <- data.frame(var = factor(names(importance), levels = names(importance), ordered = TRUE),
                   importance = importance/max(importance))
  df <- df[df$importance > remove_threshold,]

  # plot
  p <- ggplot2::ggplot(data = df, aes(x = var, y = importance)) +
    ggplot2::geom_bar(stat="identity", fill="steelblue") +
    # scale_y_continuous(trans='log10') +
    ggplot2::theme_minimal() +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Variable", y = "Importance (%)")
  p

}
