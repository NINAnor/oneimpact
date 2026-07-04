#' Predict from a bag of models or from formula + coefficients
#'
#' This function predicts responses for new data using either:
#' - a `bag` object (via `predict.bag`), or
#' - a model `formula` plus `coefs` and `weights` (via `predict.formula`).
#'
#' Predictions can be computed for full covariate data or for focal predictors only,
#' while non-focal predictors are set to a baseline values (`median`, `mean`, or `zero`).
#' It also supports ZOI-distance inputs that are internally transformed into ZOI predictors.
#'
#' @param x `[bag,list or formula]` \cr A bag of models from [oneimpact::bag_models()],
#'   or a `formula` used to build the model matrix for prediction.
#' @param newdata `[data.frame]` \cr New data for prediction. Can contain all model
#'   variables or only focal variables.
#' @param data `[data.frame=NULL]` \cr Original data used for model fitting. Used in
#'   `predict.bag()` to recover baseline/categorical levels. Not used by `predict.formula()`.
#' @param type `[character="linear"]{"linear", "exponential", "exp", "logit", "cloglog"}` \cr
#'   Prediction scale.
#' @param wmean `[logical(1)=TRUE]` \cr Whether to compute and return weighted mean prediction.
#' @param wq_probs `[numeric,vector=c(0.025, 0.5, 0.975)]` \cr Weighted quantile
#'   probabilities for prediction summaries. If `NULL`, quantiles are not returned.
#' @param include `[character="all"]` \cr Terms to include in prediction.
#'   Use `"all"` or one/more string patterns matching term names.
#' @param baseline `[character="median"]{"median", "mean", "zero"}` \cr Baseline
#'   values for non-focal predictors.
#' @param zoi `[logical(1)=FALSE]` \cr If `TRUE`, columns in `newdata` are treated as
#'   distance inputs and transformed into ZOI predictors.
#' @param zoi_shape `[character="exp_decay"]{"exp_decay", "gaussian_decay", "linear_decay", "threshold_decay"}` \cr
#'   ZOI decay shape used when `zoi = TRUE`.
#' @param which_cumulative `[character(1)="cumulative"]` \cr Pattern used to identify
#'   cumulative ZOI terms.
#' @param type_feature `[character="point"]{"point", "line", "area"}` \cr Feature type
#'   for ZOI prediction.
#' @param type_feature_recompute `[logical(1)=FALSE]` \cr Whether to recompute line-
#'   and area-feature geometry approximation for ZOI calculations.
#' @param n_features `[numeric(1)=1]` \cr Number of features used in ZOI prediction.
#' @param zoi_limit `[numeric(1)=0.05]` \cr Lower influence threshold used by non-vanishing
#'   ZOI decay functions (relevant for line-feature ZOI transformation).
#'   See [oneimpact::zoi_functions()].
#' @param resolution `[numeric(1)=100]` \cr Raster resolution used for line-feature ZOI
#'   approximation. Used when recomputing the ZOI variables for line and area features.
#' @param line_value `[numeric(1)=1]` \cr Value assigned to rasterized line cells when
#'   `type_feature = "line"`.  Used when recomputing the ZOI variables for line and
#'   area features.
#' @param coefs `[numeric vector or matrix]` \cr Coefficients used by `predict.formula()`.
#'   If matrix, rows are term names and columns are model/resample coefficients.
#' @param weights `[numeric=1]` \cr Model weights used by `predict.formula()` when
#'   combining predictions across models.
#' @param ... \cr Additional arguments.
#'
#' @return A `data.frame` (or matrix-like object) with predicted values.
#' Output columns depend on `wmean` and `wq_probs`:
#' \itemize{
#'   \item weighted quantiles (if `wq_probs` is not `NULL`)
#'   \item weighted mean (if `wmean = TRUE`)
#'   \item individual model predictions (if `wmean = FALSE` and `wq_probs = NULL`)
#' }
#'
#' @seealso [oneimpact::plot_response()], [oneimpact::create_linear_feature_zoi()].
#'
#' @example examples/bag_predict_example.R
#'
#' @export
predict <- function(x, newdata, ...) {
  UseMethod("predict")
}

#' @rdname predict
#' @export
predict.bag <- function(x,
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
                        type_feature_recompute = FALSE,
                        n_features = 1,
                        zoi_limit = 0.05, # for computing line ZOI
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
                                                                 type_feature_recompute = type_feature_recompute,
                                                                 zoi_limit = zoi_limit,
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
  pred <- predict(x$formula,
                  newdata = newdata,
                  coefs = x$coef,
                  weights = x$weights,
                  type = type,
                  wmean = wmean,
                  wq_probs = wq_probs,
                  include = include)

  pred
}

#' @rdname predict
#' @export
predict.formula <- function(x,
                            newdata,
                            coefs,
                            weights = 1,
                            type = c("linear", "exponential", "exp", "logit", "cloglog")[1],
                            wmean = TRUE,
                            wq_probs = NULL,
                            include = "all",
                            ...) {

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
                                 probs = wq_probs)))
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

predict_components <- function(x, newdata, include="all"){
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


#' Creates a line feature and gets the value of their zone of influence on locations over the feature
#'
#' The function `create_linear_feature_zoi()` computes the cumulative zone of influence (ZOI) of one single
#' linear feature so it is used to correctly create predictions and response plots for
#' linear infrastructure, considering the potential responses at multiple radii.
#'
#' This function could be extended to the nearest ZOI, if needed.
#'
#' @param radii `[numeric,vector=c(100, 250, 500, 1000, 2500, 5000, 10000)]` \cr Vector of radii for which the
#' zone of influence should be computed.
#' @param type `[character="circle"]{"circle", "Gauss", "rectangle", "exp_decay", "bartlett", "threshold"}` \cr
#' Shape of the zone of influence (ZOI), Default is `circle"`. It can assume any of the
#' possible values for the argument `type` in the function [oneimpact::dist_decay()].
#' @param radius_max `[numeric=max(radii)]` \cr Maximum radius, used to set the size of the
#' landscape/raster for ZOI computations.
#' @param res `[numeric(1)=100]` \cr Resolution for the raster created. This might impact what are the values observed
#' in the ZOI.
#' @param line_value `[numeric(1)=1]` \cr Value set to the raster line created. Default is 1.
#' It could be changed to different values if we want to represent e.g. the value in the
#' linear feature as the roads traffic or another value for spatio-temporally dynamic variables.
#'
#' @seealso [oneimpact::predict()]
#'
#' @examples
#' # create feature
#' create_linear_feature_zoi(radii = c(100, 250, 500, 1000, 2500, 5000, 10000),
#'                           type = "exp_decay",
#'                           res = 100)
#' create_linear_feature_zoi(type = "exp_decay",
#'                           type_feature_recompute = FALSE)
#'
#' @keywords internal
#' @export
create_linear_feature_zoi <- function(radii = c(100, 250, 500, 1000, 2500, 5000, 10000),
                                      type = c("circle", "Gauss", "rectangle", "exp_decay", "bartlett", "threshold",
                                               "mfilter")[1],
                                      radius_max = max(radii),
                                      zoi_limit = 0.05,
                                      type_feature_recompute = FALSE,
                                      res = 100,
                                      value = 1) {

  # compute
  if(type_feature_recompute) {
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
    # multiple the line by value
    r_line <- value * terra::rasterize(line, rr, touches= TRUE, background = 0)
    r_line
    # plot(r_line)

    # compute zois
    zois <- calc_zoi_cumulative(r_line, radius = radii, type = type, zoi_limit = zoi_limit)
    # plot(zois)

    # return(setNames(zois[radius_max/res+1, radius_max/res], radii))
    return(setNames(unlist(terra::global(zois, max)), radii))
  } else {

    # if we want to use values already computed
    radii <- c(100, 250, 500, 1000, 2500, 5000, 10000)

    if(type %in% c("bartlett", "linear", "bartlett_decay", "linear_decay")) {
      vals <- c(1, 2.6, 5, 10, 25, 50, 100)
    }
    if(type %in% c("exp_decay", "exp", "exponential_decay", "exponential")) {
      vals <- c(1.105250, 1.861974, 3.426254, 6.690854, 16.580199, 33.088243, 63.428301)
    }
    if(type %in% c("circle", "threshold")) {
      vals <- c(3, 5, 11, 21, 51, 101, 200)
    }
    if(type %in% c("Gauss")) {
      vals <- c(1.100013, 2.560138, 5.120141, 10.236728, 25.564112, 51.106633, 100.932931)
    }

    return(setNames(vals, radii))
  }


}

