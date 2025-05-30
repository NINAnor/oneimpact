#' Prediction of a bag of models to new data
#'
#' @param x `[bag,list]` \cr  A bag of models, resulting from a call to [oneimpact::bag_models()].
#' @param newdata \cr New data set to be used for prediction.
#' @param inclue vector with 1s and 0s about which variables should be included?
#'
#' @export
bag_predict <- function(x,
                        newdata,
                        type = c("linear", "exponential", "exp", "logit", "cloglog")[1],
                        wMean = T,
                        wq_probs = NULL,
                        include = "all") {

  # response
  case <- oneimpact::extract_response_strata(x$formula)$response
  # trick, check LATER ###########
  # attr(mm, "assign") |> duplicated()
  # if contrasts and NA
  # mm[is.na(mm)] <- 0

  # coefs
  coefs <- x$coef

  # subset of variables to be included
  if (include[1] != "all"){

    # terms
    include_terms <- paste0(include, collapse = "|")

    # subset coefficients
    coefs <- coefs[grepl(include_terms, rownames(coefs)), ]

    # re-set model matrix
    form_parts <- as.character(x$formula_no_strata)[3] |>
      gsub(pattern = "\n", replacement = "") |>
      strsplit(split = "+", fixed = TRUE)
    form_parts_subset <- sapply(form_parts, grep, pattern = include_terms, value = TRUE)
    form_include <- as.formula(paste0(case, " ~ -1 + ", paste(form_parts_subset, collapse = " + ")))

    nd <- newdata[, c(case, grep(include_terms, colnames(newdata), value = TRUE))]
    mm <- model.matrix(form_include, nd)
  } else {
    # model matrix
    mm <- model.matrix(x$formula_no_strata, newdata)
  }

  # prediction
  pred <- mm %*% coefs

  # if wq_probs are provided, the weighted quantiles are computed
  # if wMean is providade, the weighted mean is computed
  # if none is provided, the raw prediction is shown
  if (!is.null(wq_probs)){
    preddf <- data.frame(t(apply(pred, 1, DescTools::Quantile,
                                 weights = x$weights,
                                 type = 5,
                                 probs=wq_probs)))
    # preddf <- data.frame(t(apply(pred, 1, modi::weighted.quantile,
    #                              w = x$weights,
    #                              prob=wq_probs)))

    ########### error here
    # 50: In regularize.values(x, y, ties, missing(ties), na.rm = na.rm) :
    #   collapsing to unique 'x' values
    # when using weights = x$weights
    colnames(preddf) <- paste0("quantile:", wq_probs)
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
  if (type == "exponential" | type == "exp") {
    preddf <- exp(preddf)

  } else {
    if (type == "logit") {
      preddf <- 1 / (1 + exp(-preddf))#log(preddf/(1-preddf))
    } else {
      if (type == "cloglog") {
        preddf <- 1 - exp(-exp(preddf))#log(preddf/(1-preddf))
      }
    }
  }

  # return prediction
  return(preddf)
}

bag_predict_components <- function(x, newdata, include="all"){
  mm <- model.matrix(x$formula_no_strata, newdata)
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
