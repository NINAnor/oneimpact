#' Prediction of a bag of models to new data
#'
#' The function `bag_predict` makes a prediction from
#'
#' @param x `[bag,list or formula]` \cr  A bag of models, resulting from a call to [oneimpact::bag_models()],
#' or a `formula` used to fit the models in the bag.
#' @param newdata \cr New data set to be used for prediction. It can include all the variables in the formula
#' or only those for which the user is interested in making a prediction from.
#' @param inclue vector with 1s and 0s about which variables should be included?
#'
#' @example examples/bag_predict_examople.R
#'
#' @export
bag_predict <- function(x,
                        newdata,
                        type = c("linear", "exponential", "exp", "logit", "cloglog")[1],
                        wmean = TRUE,
                        wq_probs = NULL,
                        include = "all", ...) {
  UseMethod("bag_predict")
}

#' @rdname bag_predict
#' @export
bag_predict.bag <- function(x,
                            newdata,
                            data = NULL,
                            type = c("linear", "exponential", "exp", "logit", "cloglog")[1],
                            wmean = TRUE,
                            wq_probs = NULL,
                            include = "all",
                            baseline = c("median", "mean", "zero")[1],
                            zoi = FALSE,
                            zoi_shape = c("exp_decay", "gaussian_decay", "linear_decay", "threshold_decay")[1],
                            which_cumulative = "cumulative",
                            type_feature = c("point", "line", "area")[1],
                            n_features = 1,
                            resolution = 100, # resolution for the raster created in create_line_feature_zoi
                            line_value = 1, # value sey to the linear inftastructure raster created in create_line_feature_zoi
                            ...) {

  # store new data entered
  dfvar <- newdata

  # baselines for plotting and predicting!! mean? median?
  # define the baseline
  bs <- baseline
  if (baseline[1] == "median"){
    baseline <- x$data_summary[rownames(x$data_summary) == "50%",]
  } else {
    if(baseline[1] == "mean") {
      baseline <- x$data_summary[rownames(x$data_summary) == "mean",]
    } else {
      if (baseline[1] == "zero"){
        baseline <- x$data_summary[5,]
        baseline[1, 1+which(x$numeric_covs)] <- 0 # zero for numeric ones
        if(!zoi) baseline[1, 1+which(!x$numeric_covs)] <- sapply(unname(which(!x$numeric_covs)),
                                                                 function(z) unique(data[,names(x$numeric_covs[z])])[1])
      } else {
        stop(paste0("Invalid 'baseline' parameter: ", baseline, ". Pleas re-define."))
      }
    }
  }

  # is the variable a ZOI variable?
  if (!zoi){

    # new data
    newdata <- baseline
    newdata <- newdata[rep(1, nrow(dfvar)),]
    newdata[, names(dfvar)] <- dfvar

  } else{

    # ZOI radii
    zois <- names(x$data_summary)[grep(names(dfvar)[1], names(x$data_summary))]
    # zois <- as.numeric(gsub(names(dfvar)[1], "", zois))
    zoi_radii <- as.numeric(gsub("\\D", "", zois))

    # compute ZOI for the intended distances
    dfvar2 <- dfvar[, rep(1, length(zoi_radii)), drop = FALSE]
    is_cumulative <- grepl(pattern = which_cumulative, zois)
    fact <- rep(1, length(zoi_radii))
    if(type_feature == "line") fact <- create_linear_feature_zoi(radii = zoi_radii,
                                                                 type = zoi_shape,
                                                                 radius_max = max(zoi_radii),
                                                                 res = resolution,
                                                                 value = line_value)
    dfvar2 <- as.data.frame(do.call("cbind", lapply(c(1:ncol(dfvar2)), function(i) {
      n_feat <- ifelse(is_cumulative[i], n_features, 1)
      n_feat*fact[i]*oneimpact::dist_decay(dfvar2[,i], radius = zoi_radii[i], type = zoi_shape) })))
    names(dfvar2) <- zois#paste0(names(dfvar)[1], zoi_radii)

    # cumulative vars
    #dfvar2[,grepl("cumulative", colnames(dfvar2))]*100

    # new data
    newdata <- baseline
    newdata <- newdata[rep(1, nrow(dfvar)),]
    newdata[, names(dfvar2)] <- dfvar2
    #if (ncol(dfvar)==2){newdata[,names(dfvar)[2]] <- dfvar[,2]}
  }

  # make sure prediction works even if categorical variables are constant
  # get that from the summary instead of dat, where from? maybe a new object
  if(!(zoi & bs == "zero")) { # other variables are ignored in the case of zoi plots
    for(i in which(!x$numeric_covs)) {
      if(class(data[,names(x$numeric_covs[i])]) == "factor") {
        newdata[,i+1] <- factor(newdata[,i+1], levels = levels(data[,names(x$numeric_covs[i])])) # factor
      } else {
        newdata[,i+1] <- factor(newdata[,i+1], levels = sort(unique(data[,names(x$numeric_covs[i])]))) # character
      }
    }
  }

  # predict for new data set

  # predict only for those variables/coefs of interest for ZOI variables
  include <- if(zoi) colnames(dfvar) else include
  # predict
  pred <- bag_predict(x$formula,
                      newdata = newdata,
                      coefs = x$coef,
                      weights = x$weights,
                      type = type,
                      wmean = wmean,
                      wq_probs = wq_probs,
                      include = include)

  pred
}


#' @param x `[formula]` \cr Formula used to fit the model.
#' @param coefs `[vector,numeric]` \cr Either a named vector of coefficients (in case there is only one
#' model) or a matrix of coefficients, with rownames as the term names and columns as the different models/resamples.
#' @param weights `[vector,numeric=1]` \cr Vector of weights for the different models/resamples, i.e.
#' the column from the `coefs` object with coefficients. A single number (by default, 1) in case there is only
#' one model (`coefs` is a vector).
#'
#' @rdname bag_predict
#' @export
bag_predict.formula <- function(x,
                                newdata,
                                coefs,
                                weights = 1,
                                type = c("linear", "exponential", "exp", "logit", "cloglog")[1],
                                wmean = TRUE,
                                wq_probs = NULL,
                                include = "all") {

  # formula with not strata
  wcols <- extract_response_strata(x, covars = TRUE)
  formula_no_strata <- as.formula(paste0(wcols$response, " ~ -1+", wcols$covars))

  # repeat weights == 1 if necessary
  if(length(weights) == 1 & is.matrix(coefs)) {
    weights <- coefs[1,]
    weights[] <- 1/length(weights)
  }

  # coefs
  # coefs <- x$coef

  # subset of variables to be included
  if (include[1] != "all"){

    # terms
    include_terms <- paste0(include, collapse = "|")
    # subset coefficients
    if(is.matrix(coefs)) {
      coef_names <- rownames(coefs)
      coefs <- coefs[grepl(include_terms, coef_names), ]
    } else {
      coef_names <- names(coefs)
      coefs <- coefs[grepl(include_terms, coef_names)]
    }

    # re-set model matrix
    form_parts <- as.character(formula_no_strata)[3] |>
      gsub(pattern = "\n", replacement = "") |>
      strsplit(split = "+", fixed = TRUE)
    form_parts_subset <- sapply(form_parts, grep, pattern = include_terms, value = TRUE)
    form_include <- as.formula(paste0(" ~ -1 + ", paste(form_parts_subset, collapse = " + ")))

    nd <- newdata[, c(grep(include_terms, colnames(newdata), value = TRUE)), drop = FALSE]
    mm <- model.matrix(form_include, nd)
  } else {
    # model matrix
    mm <- model.matrix(formula_no_strata, newdata)
  }

  # prediction
  pred <- mm %*% coefs

  # if wq_probs are provided, the weighted quantiles are computed
  # if wmean is providade, the weighted mean is computed
  # if none is provided, the raw prediction is shown
  if (!is.null(wq_probs)){
    preddf <- data.frame(t(apply(pred, 1, DescTools::Quantile,
                                 weights = weights,
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
    if (wmean){ preddf$mean <- as.vector(pred %*% weights) }
  }else{
    if (wmean){
      preddf <- data.frame(mean = as.vector(pred %*% weights))
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


#' The function `create_linear_feature_zoi` created ZOI of one single linear feature so it is used to
#' correctly create response plots for linear infrastructure, considering the potential
#' responses at multiple radii.
#'
#' @export
create_linear_feature_zoi <- function(radii = c(100, 250, 500, 1000, 2500, 5000, 10000),
                                      type = c("circle", "Gauss", "rectangle", "exp_decay", "bartlett", "threshold",
                                               "mfilter")[1],
                                      radius_max = 10000,
                                      res = 100,
                                      value = 1) {

  # create line feature
  line <- data.frame(
    x = c(-radius_max, radius_max),
    y = c(0, 0)
  ) |>
    terra::vect(geom = c("x", "y"), keepgeom = TRUE) |>
    terra::as.lines()
  # plot(line)
  # rasterize it
  rr <- terra::rast(nrows = radius_max, ncols = radius_max,
                    xmin = -radius_max, xmax = radius_max,
                    ymin = -radius_max, ymax = radius_max,
                    res = res)
  r_line <- value * terra::rasterize(line, rr, touches= TRUE, background = 0)
  r_line
  # plot(r_line)

  # compute zois
  zois <- calc_zoi_cumulative(r_line, radius = radii, type = type)
  # plot(zois)

  # return(setNames(zois[radius_max/res+1, radius_max/res], radii))
  return(setNames(unlist(global(zois, max)), radii))
}
