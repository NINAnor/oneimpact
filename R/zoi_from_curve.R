#' Get estimates of zone of influence (ZOI) from response curves
#'
#' This generic function computes ZOI metrics (maximum effect size,
#' ZOI radius and impact) for ZOI predictor variables based on response curves.
#' The ZOI radius is estimated as the distance/radius in which the
#' relative selection strength corresponds to on a given percentage of
#' the maximum effect size (e.g. 95% ZOI radius). The impact accounts for
#' both the effect size and the ZOI radius and corresponds to the area under
#' (or over, if negative) the ZOI response curve.
#' The function supports two types of input: a data.frame of predictions
#' or a bag of models.
#'
#' @param x Either a `data.frame` containing response curve predictions for a single variable,
#'   or a `bag` object containing an ensemble of models.
#' @param ... Additional arguments passed to the appropriate method.
#'
#' @return A `data.frame` or a `list` containing ZOI metrics:
#' - `max_effect_size`: Maximum effect size on the relative selection strength (y axis).
#' - `zoi_radius`: Distance at which the effect drops below a threshold,
#' defined by the parameter `percentage`.
#' - `effect_zoi_radius`: Relative selection strength (y axis) valye in which the ZOI
#' is reached.
#' - `impact`: Area under the curve up to the ZOI radius,
#' combining the varying effect size with distance.
#' Each ZOI measure presents mean, median, CI lower, and CI upper.
#'
#' @seealso [oneimpact::predict()], [oneimpact::plot_response()], [oneimpact::weirdness()]
#' ##example examples/zoi_from_curve_example.R
#'
#' @name zoi_from_curve
#' @export
zoi_from_curve <- function(x, ...) {
  UseMethod("zoi_from_curve")
}

#' @param percentage `[numeric(1)=0.95]` \cr Numeric between 0 and 1. Defines the
#' threshold for ZOI radius as a proportion of the maximum effect size.
#' Default is `0.95`.
#' @param curve `[character(1)=c("mean", "median")]` \cr Character vector.
#' Which central tendency curves to use: `"median"`, `"mean"`, or both.
#' @param ci `[logical(1)=TRUE]` \cr Logical. Whether to compute ZOI estimates for
#' the upper and lower limits of the confidence interval. Default is TRUE.
#' @param type `[character(1)="linear"]{"linear", "exp"}` \cr Character. Defines whether
#' the calculation of ZOI should be based on the prediction of at linear
#' or response (exponential) scale: `"linear"` or `"exp"`, respectively.
#' @param mean_col_name `[character="mean"]` \cr Name of the column containing
#' the mean response curve.
#' @param median_col_name `[character="quantile:0.5"]` \cr Name of the column
#' containing the median response curve.
#' @param ci_col_name `[character=c("quantile:0.255", "quantile:0.975")]` \cr Character
#' vector of length 2. Names of columns for lower and upper confidence intervals.
#'
#' @rdname zoi_from_curve
#' @export
zoi_from_curve.data.frame <- function(x,
                                      percentage = 0.95,
                                      curve = c("median", "mean"),
                                      ci = TRUE,
                                      type = c("linear", "exp")[1],
                                      # n_features = 1,
                                      mean_col_name = "mean",
                                      median_col_name = "quantile:0.5",
                                      ci_col_name = c("quantile:0.025", "quantile:0.975")) {

  # get predictor / ZOI variable
  xvar <- colnames(x)[1]

  if(type == "linear") {
    ref <- 0
  } else {
    ref <- 1
  }

  # initialize output
  # mean, median, 0.025, 0.975
  max_effect_size <- rep(NA, 4)
  zoi_radius <- rep(NA, 4)
  effect_zoi_radius <- rep(NA, 4)
  impact <- rep(NA, 4)

  # main and median effects
  i <- 1
  for(i in seq_along(curve)) {

    # get values of the main response - either mean or median
    if(curve[i] == "median") {
      id <- 2
      main_response <- x[[median_col_name]]
    } else {
      id <- 1
      main_response <- x[[mean_col_name]]
    }

    # # max effect size
    # if(type == "linear") {
    # max_effect_size_main <- main_response[which.max(abs(main_response))]
    # } else {
    max_effect_size_main <- main_response[which.max(abs(main_response - ref))]
    # }

    # zoi radius main
    # y_value_percentage <- (1 - percentage) * (max_effect_size_main)
    # x_radius_index <- which(abs(main_response) < abs(y_value_percentage))[1] - 1
    # zoi_radius_main <- x[[xvar]][x_radius_index]

    y_value_percentage <- (1 - percentage) * (max_effect_size_main - ref)
    if(type == "linear") {
      x_radius_index <- which(abs(main_response) < abs(y_value_percentage + ref))[1] - 1
    } else {
      x_radius_index <- which(abs(main_response) > abs(y_value_percentage + ref))[1] - 1
    }
    if(is.na(x_radius_index)) {
      x_radius_index <- length(main_response)
    }
    zoi_radius_main <- x[[xvar]][x_radius_index]

    # impact
    x_vals <- x[[xvar]][1:x_radius_index]
    y_vals <- main_response[1:x_radius_index] - ref
    signal <- ifelse(y_vals[1] < 0, -1, 1)
    y_vals <- abs(y_vals) - min(abs(y_vals)) # get positive and discount are above y(ZOI)

    impact_main <- signal * DescTools::AUC(x_vals, y_vals)

    max_effect_size[id] <- max_effect_size_main
    zoi_radius[id] <- zoi_radius_main
    effect_zoi_radius[id] <- y_value_percentage + ref
    impact[id] <- impact_main
  }

  # repeat that for the CI
  if(ci) {

    ci_id <- 2
    for(ci_id in seq(ci_col_name)) {

      ci_response <- x[[ci_col_name[ci_id]]]

      # max effect size
      ci_max <- ci_response[which.max(abs(main_response - ref))]

      # zoi radius main
      # y_value_percentage <- (1 - percentage) * max_effect_size
      if(type == "linear") {
        x_radius_ci_index <- which(abs(ci_response) < abs(y_value_percentage + ref))[1] - 1
      } else {
        x_radius_ci_index <- which(abs(ci_response) > abs(y_value_percentage + ref))[1] - 1
      }
      if(is.na(x_radius_ci_index)) {
        x_radius_ci_index <- length(main_response)
      }
      zoi_radius_ci <- x[[xvar]][x_radius_ci_index]

      if(length(zoi_radius_ci) == 0) {
        zoi_radius_ci <- NA
      }

      # impact
      x_vals <- x[[xvar]][1:x_radius_ci_index]
      y_vals <- ci_response[1:x_radius_ci_index] - ref
      signal <- ifelse(y_vals[1] < 0, -1, 1)
      y_vals <- abs(y_vals) - min(abs(y_vals)) # get positive and discount are above y(ZOI)

      impact_ci <- signal * DescTools::AUC(x_vals, y_vals)

      max_effect_size[ci_id + 2] <- ci_max
      zoi_radius[ci_id + 2] <- zoi_radius_ci
      effect_zoi_radius[ci_id + 2] <- y_value_percentage + ref
      impact[ci_id + 2] <- impact_ci
    }
  }

  out <- data.frame(max_effect_size = max_effect_size,
              zoi_radius = zoi_radius,
              effect_zoi_radius = effect_zoi_radius,
              impact = impact) |>
    # dplyr::bind_rows() |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "zoi_measure")

  colnames(out)[-1] <- c(mean_col_name, median_col_name, ci_col_name[1], ci_col_name[2])
  out
}

#' @param data `[data.frame]` \cr The original dataset used for model fitting.
#' @param include `[character="all"]` \cr Character. Either `"all"` or a
#' regex pattern to filter selected ZOI variables.
#' @param return_predictions `[logical=FALSE]` \cr Logical. Whether to return
#' the prediction curves alongside ZOI metrics. If `TRUE`, the output is necessarily
#' a `list` with predictions and the ZOI parameters.
#' @param return_format `[character="df"]{"list", "df"}` \cr
#' Format of the returned ZOI metrics. Either a list of data.frames (if `return_format = "list"`),
#' one for each variable, or a single `data.frame` (default, if `return_format = "df"`).
#' @param wq_probs `[numeric,vector=c(0.025, 0.975)]` \cr Numeric vector of quantiles
#' used for prediction summaries.
#' @param n_features `[numeric=1]` \cr Number of features used in ZOI prediction.
#' It can a single number (considered the same for all ZOI variables) or a vector
#' with the same number of elements as ZOI variables in the model.
#' @param radius_max `[numeric=NULL]` \cr Numeric. Maximum distance/radius to use for
#' prediction curves. If `NULL` (default), the maximum value present in the bag's
#' predictor table is used.
#' @param baseline `[character="zero"]` \cr Character. Baseline used in `predict()` (e.g., `"zero"`).
#' @param type_feature `[character="point"]` \cr Character or vector. Type of spatial feature used in
#' `predict()`.
#' @param type_feature_recompute `[logical=FALSE]` \cr Logical. Whether to recompute spatial
#' features within `predict()`, for linear features.
#' @param resolution `[numeric=200]` \cr Integer. Resolution used in the recomuptation
#' of ZOIs for linear features.
#' @param radii `[vector]` \cr Numeric vector. Radii used for ZOI modeling.
#' @param zoi_shape `[character]` \cr Character. Shape of the ZOI used in the model
#' (e.g., `"circle"`, `"Gauss"`, `"exp_decay"`).
#'
#' @return If `x` is a bag object, the function returns wither a `list` or
#' `data.frame` of ZOI measures for each ZOI variable in the bag.
#' If `return_predictions = TRUE`, also returns the prediction curves.
#'
#' @rdname zoi_from_curve
#' @export
zoi_from_curve.bag <- function(x,
                               data,
                               include = "all",
                               percentage = 0.95,
                               curve = c("median", "mean"),
                               type = c("linear", "exp")[1],
                               return_predictions = FALSE,
                               return_format = c("list", "df")[2],
                               ci = TRUE,
                               wq_probs = c(0.025, 0.5, 0.975),
                               n_features = 1,
                               mean_col_name = "mean",
                               median_col_name = "quantile:0.5",
                               ci_col_name = c("quantile:0.025", "quantile:0.975"),
                               radius_max = NULL,
                               baseline = "zero",
                               type_feature = "line",
                               type_feature_recompute = TRUE,
                               resolution = 200,
                               radii = c(100, 250, 500, 1000, 2500, 5000, 10000),
                               zoi_shape = c("circle", "Gauss", "rectangle", "exp_decay", "bartlett", "threshold",
                                             "mfilter")[1],
                               ...) {

  # get ZOI variables and terms from predictor table
  pred_table <- x$parms$predictor_table
  zoi_vars <- pred_table$variable[pred_table$is_zoi == 1]
  zoi_terms <- pred_table$term_zoi[pred_table$is_zoi == 1]
  zoi_radii <- pred_table$zoi_radius[pred_table$is_zoi == 1]

  # set radius from predictor table, if null
  if(is.null(radius_max)) radius_max <- max(zoi_radii)

  # unique variables
  zoi_vars_unique <- unique(zoi_vars)
  if(include != "all") {
    vars <- grep(include, zoi_vars_unique, value = TRUE)
    if(length(vars) < 1) {
      stop(paste0("Variable(s) ", paste0(include, collapse = ","), " not present in the bag."))
    }
    zoi_vars_unique <- vars
  }

  # check parameters
  if(length(type_feature) == 1) {
    type_feature <- rep(type_feature, times = length(zoi_vars_unique))
  }

  if(length(n_features) == 1) {
    n_features <- rep(n_features, times = length(zoi_vars_unique))
  }

  # compute predictions
  i <- 3
  dfs <- lapply(seq_along(zoi_vars_unique), function(i) {
    dfvar <- data.frame(var = seq(0, radius_max, length.out = 10001))
    names(dfvar) <- zoi_vars_unique[i]
    type_feat <- type_feature[i]
    zoi_r <- zoi_radii[zoi_vars == zoi_vars_unique[i]]
    pred <- oneimpact::predict(x,
                               newdata = dfvar,
                               data = data,
                               type = type,
                               wq_probs = wq_probs,
                               zoi = TRUE,
                               n_features = n_features[i],
                               baseline = baseline,
                               type_feature = type_feat,
                               type_feature_recompute = type_feature_recompute,
                               resolution = resolution,
                               radii = zoi_r,
                               zoi_shape = zoi_shape, ...)
    cbind(dfvar, pred)
  })
  names(dfs) <- zoi_vars_unique

  # compute zoi
  i <- 1
  zois <- lapply(seq_along(dfs), function(i) {
    zoi_from_curve(dfs[[i]], type = type,
                   percentage = percentage,
                   curve = curve,
                   ci = ci,
                   mean_col_name = mean_col_name,
                   median_col_name = median_col_name,
                   ci_col_name = ci_col_name)
  })
  names(zois) <- zoi_vars_unique

  if(return_format == "df") {
    # i <- 1
    zois <- lapply(seq_along(zois), function(i) {
      zois[[i]] |>
        dplyr::mutate(variable = zoi_vars_unique[i])
    }) |>
      dplyr::bind_rows() |>
      dplyr::relocate(variable, .before = 1) |>
      tibble::as_tibble()
  }

  # return tables with zois
  if(return_predictions) {
    return(list(predictions = df, zoi = zois))
  } else {
    return(zois)
  }
}
# implement function var
# implement function bag - all vars


