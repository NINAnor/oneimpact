#' Prediction of a bag of models to new data
#'
#' @param x `[list]` \cr  A bag of models, resulting from a call to [oneimpact::bag_models()].
#' @param newdata \cr New data set to be used for prediction.
#' @param inclue vector with 1s and 0s about which variables should be included?
#'
#' @export
bag_predict <- function(x,
                        newdata,
                        type = c("linear", "exponential")[2],
                        wMean = T,
                        wQ_probs = NULL,
                        include = "all") {

  # model matrix
  mm <- model.matrix(x$mm_formula, newdata)
  # trick, check LATER ###########
  # attr(mm, "assign") |> duplicated()
  # if contrasts and NA
  # mm[is.na(mm)] <- 0

  # coefs
  coefs <- x$coef

  # subset of variables to be included
  if (include[1]!="all"){
    coefs <- coefs * include
  }

  # prediction
  pred <- mm %*% coefs

  # if wQ_probs are provided, the weighted quantiles are computed
  # if wMean is providade, the weighted mean is computed
  # if none is provided, the raw prediction is shown
  if (!is.null(wQ_probs)){
    preddf <- data.frame(t(apply(pred, 1, DescTools::Quantile, #weights = NULL,
                                 weights = x$weights,
                                 type = 5,
                                 probs=wQ_probs)))
    ########### error here
    # 50: In regularize.values(x, y, ties, missing(ties), na.rm = na.rm) :
    #   collapsing to unique 'x' values
    # when using weights = x$weights
    colnames(preddf) <- paste0("quantile:", wQ_probs)
    rownames(preddf) <- NULL
    if (wMean){ preddf$mean <- as.vector(pred %*% x$weights) }
  }else{
    if (wMean){
      preddf <- data.frame(mean = as.vector(pred %*% x$weights))
    }else{
      preddf <- pred
    }
  }

  # should result be in linear or exp scale?
  if (type=="exponential") { preddf <- exp(preddf) }

  # return prediction
  return(preddf)
}

bag_predict_components <- function(x, newdata, include="all"){
  mm <- model.matrix(x$mm_formula, newdata)
  coefs <- x$coef %*% x$weights
  if (include[1]!="all"){
    coefs <- coefs * include
  }
  wrow <- which(coefs!=0)
  pred <- as.data.frame(lapply(wrow, function(x, mm, coefs){mm[,x] * coefs[x]}, mm=mm, coefs=coefs))
  names(pred) <- x$var_names[wrow]
  return(pred)
}

combine_zoi_components <- function(x, zoi=c(250, 500, 1000, 2500, 5000, 10000)){
  zoi <- zoi[order(-zoi)]
  zoi <- paste0("_", zoi)
  for (i in zoi){names(x) <- gsub(i, "", names(x))}

  zoi_vars <- unique(names(x)[duplicated(names(x))])
  y <- x[,!duplicated(names(x))]
  for (i in zoi_vars){y[,names(y)==i] <- rowSums(x[,names(x)==i])}
  return(y)
}
