#' Computes ecological weirdness for a fitted model or it estimated  coefficients
#'
#' @param x Bag.
#' @param data `[data.frame]` \cr The original, complete data used for model fitting.
#' @param measure `[string(1)]{""coef_sign", "n_crosses", "response_area""}` \cr Measure used
#' to quantify "weirdness" in the model or coefficients, based on the coeffcients and the response
#' plots for each type of covariate with zone of influence in a model.
#' It can be one or multiple of these options:
#' - `"coef_sign"`: counts the number of coefficients whose sign is opposite to the ecologically expected sign;
#' - `"n_crosses"`: counting the number of crosses in sign for the coefficients of the same covariate;
#' - `"response_area"`: computing the area under the response plot curve which is in the unexpected
#' direction.
#' @param which_coef \cr Which measure to use for the coefficients, when `measure = "coef_sign"`. If `count` (default),
#' only the sign matterns and we count the number of coefficients with unexpected sign.
#' If `sum`, we count the sum of the (standardized) coefficients, to also account for their magnitude.
#' @param expected_sign `[numeric(1)=-1]` \cr Expected sign of the coefficient. Either -1 (negative),
#' +1 (positive), or 0 (no effect).
#' @param zero_coefficient_limit `[numeric(1)=1e8]` \cr Value above which an estimated coefficient is considered
#' non-zero. Default is 1e-8. Useful for comparing coefficients which are expected to be zero (i.e. to have no effect).
#'
#' @example examples/weirdness_example.R
#'
#' @name weirdness
#' @export
weirdness <- function(x, ...) {
  UseMethod("weirdness")
}

# here x is a numeric vector of coefficients
#' @rdname weirdness
#' @export
weirdness.numeric <- function(x,
                              which_coef_sign = c("count", "sum", "raw", "index")[1],
                              expected_sign = -1,
                              zero_coefficient_limit = 1e-8) {

  # x = vector of standardized coefficients

  # for count of unexpected signs
  if(which_coef_sign == "count") {
    weird <- ifelse(expected_sign == 0,
                    sum(abs(x) > zero_coefficient_limit),
                    sum(x*expected_sign < 0))
  } else {

    if(which_coef_sign == "sum") {
      # for sum of unexpected signs
      weird <- ifelse(expected_sign == 0,
                      sum(abs(x[abs(x) > zero_coefficient_limit])),
                      abs(sum(x[x*expected_sign < 0])))
    } else {

      if(which_coef_sign == "raw") {
        # to identify each term
        weird <- if(expected_sign == 0)
          x[which(abs(x) > zero_coefficient_limit)] else
            x[which(x*expected_sign < 0, )]
      } else {

        if(which_coef_sign == "index") {
          # to identify each term
          weird <- if(expected_sign == 0)
            which(abs(x) > zero_coefficient_limit) else
              which(x*expected_sign < 0, )
        } else {
          stop("The argument 'which_coef_sign' is invalid, please check.")
        }

      }
    }

  }

  # return
  weird
}

#' @rdname weirdness
#' @export
weirdness.data.frame <- function(x,
                                 expected_sign = -1,
                                 response = c("mean", "mid")[1],
                                 measure = c("n_crosses", "where_crosses",
                                             "response_area_opposite", "response_area_ratio",
                                             "n_inflection", "difference_inflection", "response_area_inflection")[1]) {

  # x = data.frame with first column as the explanatory variable and columns "mid" or "mean" as response
  resp_var <- x[[response]]

  # compute number of crosses
  if(any(grepl("n_crosses", measure))) {
    weird <- sum(sapply(2:length(resp_var), function(i) resp_var[i]/resp_var[i-1]) < 0)
  } else {

    if("where_crosses" %in% measure) {
      weird <- x[[1]][which(sapply(2:length(resp_var), function(i) resp_var[i]/resp_var[i-1]) < 0)]
      return(weird)
    }
    if(any(grepl("where_crosses_index", measure))) {
      weird <- which(sapply(2:length(resp_var), function(i) resp_var[i]/resp_var[i-1]) < 0)
      return(weird)
    }

    if(any(grepl("response_area_opposite|response_area_ratio", measure))) {
      # AUC for places on the side of the curve opposite to expectation
      area_opposite <- DescTools::AUC(x[[1]][resp_var*expected_sign < 0], resp_var[resp_var*expected_sign < 0])
      # AUC for places on the side of the curve as expected
      area_expected <- DescTools::AUC(x[[1]][resp_var*expected_sign >= 0], resp_var[resp_var*expected_sign >= 0])

      if(any(grepl("response_area_opposite", measure))) {
        weird <- area_opposite
      } else {
        weird <- area_opposite/abs(area_expected)
      }

    } else {

      if(any(grepl("n_inflection", measure))) {
        weird <- sum(inflection(resp_var))
      } else {

        if(any(grepl("difference_inflection", measure))) {
          which_inflection <- which(inflection(resp_var))
          if(length(which_inflection) < 2) {
            weird <- 0
          } else {
            weird <- sum(abs(diff(resp_var[which_inflection])))
          }

        }

        if(any(grepl("response_area_inflection", measure))) {
          stop("Response area between inflection points to be implemented.")
        }

      }

    }
  }

  # return
  ifelse(is.na(weird), 0, weird)
}

#' @export
weirdness.bag <- function(x,
                          data,
                          measure = c("coef_sign",
                                      "n_crosses", "where_crosses",
                                      "response_area_opposite",
                                      "n_inflection", "difference_inflection"),
                          wmean = TRUE,
                          which_coef_sign = c("count", "sum")[1],
                          expected_sign = -1,
                          zero_coefficient_limit = 1e-8,
                          which_n_cross = c("mean", "sum")[1],
                          response = c("mean", "mid")[1],
                          baseline = "zero",
                          type_feature_recompute = TRUE,
                          resolution = 200,
                          type_feature =  "point",
                          radii = c(100, 250, 500, 1000, 2500, 5000, 10000),
                          type = c("circle", "Gauss", "rectangle", "exp_decay", "bartlett", "threshold",
                                   "mfilter")[1],
                          radius_max = NULL,
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

  # check parameters
  if(length(type_feature) == 1) {
    type_feature <- rep(type_feature, times = length(zoi_vars_unique))
  }

  # coefs
  if(wmean) {
    coefs <- x$coef %*% x$weights
  } else {
    coefs <- x$coef[, x$weights > 0, drop = FALSE]
  }
  # coefs_std <- rescale_coefficients()

  # for now, subset the coefs for zoi variables
  coefs <- coefs[unlist(sapply(paste0(zoi_terms, "$"), function(i) grep(i, rownames(coefs)))),, drop = FALSE]

  weirdness_measures <- list(n_coefs = nrow(coefs),
                             n_resamples = x$n_above_threshold,
                             coef_sign_index = NULL,
                             coef_sign_names = NULL,
                             coef_sign_radii = NULL,
                             coef_sign_value = NULL,
                             coef_sign = NULL,
                             coef_sign_sum = NULL,
                             cross_index = NULL,
                             where_crosses = NULL,
                             n_crosses = NULL,
                             n_crosses_total = NULL,
                             response_area_opposite = NULL,
                             response_area_opposite_total = NULL,
                             response_area_ratio = NULL,
                             response_area_ratio_total = NULL,
                             n_inflection = NULL,
                             n_inflection_total = NULL,
                             difference_inflection = NULL,
                             difference_inflection_total = NULL)

  # compute measure for the sign of coefficients
  if("coef_sign" %in% measure) {
    # if(wmean) {
    # weirdness_measures$coef_sign_index <- sapply(zoi_vars_unique, function(i) {
    #   weirdness(coefs[grepl(i, rownames(coefs)),,drop = FALSE],
    #             which_coef = "index",
    #             expected_sign = expected_sign,
    #             zero_coefficient_limit = zero_coefficient_limit)
    # }, simplify = FALSE, USE.NAMES = TRUE)
    #
    # weirdness_measures$coef_sign_names <- sapply(zoi_vars_unique, function(i) {
    #   zoi_terms[grepl(i, zoi_terms)][weirdness_measures$coef_sign_index[[i]]]
    # }, simplify = FALSE, USE.NAMES = TRUE)
    #
    # weirdness_measures$coef_sign_radii <- sapply(zoi_vars_unique, function(i) {
    #   zoi_radii[grepl(i, zoi_terms)][weirdness_measures$coef_sign_index[[i]]]
    # }, simplify = FALSE, USE.NAMES = TRUE)
    #
    # weirdness_measures$coef_sign_value <- sapply(zoi_vars_unique, function(i) {
    #   coefs[grepl(i, rownames(coefs)),,drop = FALSE][weirdness_measures$coef_sign_index[[i]]]
    # }, simplify = FALSE, USE.NAMES = TRUE)
    #
    # weirdness_measures$coef_sign <- sapply(weirdness_measures$coef_sign_index, length)
    #
    # weirdness_measures$coef_sign_sum <- sum(weirdness_measures$coef_sign)
    # }

    # this code now works for only the mean (1 single model) and for the
    # bag of models
    weirdness_measures$coef_sign_index <- sapply(zoi_vars_unique, function(i) {
      lapply(seq_len(ncol(coefs)), function(z) {
        weirdness(coefs[grepl(i, rownames(coefs)), z, drop = FALSE],
                  which_coef = "index",
                  expected_sign = expected_sign,
                  zero_coefficient_limit = zero_coefficient_limit)
      })
    }, simplify = FALSE, USE.NAMES = TRUE)

    weirdness_measures$coef_sign_names <- sapply(zoi_vars_unique, function(i) {
      lapply(seq_len(ncol(coefs)), function(z) {
        zoi_terms[grepl(i, zoi_terms)][weirdness_measures$coef_sign_index[[i]][[z]]]
      })
    }, simplify = FALSE, USE.NAMES = TRUE)

    weirdness_measures$coef_sign_radii <- sapply(zoi_vars_unique, function(i) {
      lapply(seq_len(ncol(coefs)), function(z) {
        zoi_radii[grepl(i, zoi_terms)][weirdness_measures$coef_sign_index[[i]][[z]]]
      })
    }, simplify = FALSE, USE.NAMES = TRUE)

    weirdness_measures$coef_sign_value <- sapply(zoi_vars_unique, function(i) {
      lapply(seq_len(ncol(coefs)), function(z) {
        coefs[grepl(i, rownames(coefs)),,drop = FALSE][weirdness_measures$coef_sign_index[[i]][[z]]]
      })
    }, simplify = FALSE, USE.NAMES = TRUE)

    weirdness_measures$coef_sign <- sapply(weirdness_measures$coef_sign_index, function(i) sapply(i, length))

    weirdness_measures$coef_sign_sum <- sum(weirdness_measures$coef_sign)

  }
  # possibility: sum of absolute values of standardized coefficients that are against the expected sign

  # Computing data.frames with predictions for analyzing response plot curves

  # if we want to look into the wmean or wmedian curve
  if(wmean) {
    dfs <- lapply(seq_along(zoi_vars_unique), function(i) {
      dfvar <- data.frame(var = seq(0, radius_max, length.out = 10001))
      names(dfvar) <- zoi_vars_unique[i]
      type_feat <- type_feature[i]
      plot_response(x,
                    dfvar = dfvar,
                    data = data,
                    type = "linear",
                    zoi = TRUE,
                    type_feature_recompute = type_feature_recompute,
                    resolution = resolution,
                    type_feature = type_feat,
                    # resolution = 300,
                    baseline = baseline,
                    ci = TRUE,
                    indiv_pred = FALSE,
                    ggplot = FALSE,
                    ...)
    })
  } else {

    # if not, look into each individual plots
    dfs <- lapply(seq_along(zoi_vars_unique), function(i) {
      dfvar <- data.frame(var = seq(0, radius_max, length.out = 10001))
      names(dfvar) <- zoi_vars_unique[i]
      type_feat <- type_feature[i]
      plot_response(x,
                    dfvar = dfvar,
                    data = data,
                    type = "linear",
                    zoi = TRUE,
                    type_feature_recompute = type_feature_recompute,
                    resolution = resolution,
                    type_feature = type_feat,
                    # resolution = 300,
                    baseline = baseline,
                    wq_probs = NULL,
                    ci = FALSE,
                    indiv_pred = TRUE,
                    ggplot = FALSE,
                    ...) |>
        dplyr::select(-mean)
    })

  }

  names(dfs) <- zoi_vars_unique
  # str(dfs)

  # computes n crosses
  if("n_crosses" %in% measure) {

    # if(wmean) {
    #   weirdness_measures$n_crosses <- sapply(dfs, function(i) {
    #     weirdness(i,
    #               expected_sign = expected_sign,
    #               response = response,
    #               measure = "where_crosses")
    #   })
    #
    #   weirdness_measures$n_crosses <- sapply(dfs, function(i) {
    #     weirdness(i,
    #               expected_sign = expected_sign,
    #               response = response,
    #               measure = "n_crosses")
    #   })
    #   weirdness_measures$n_crosses_total <- sum(weirdness_measures$n_crosses)
    # } else {
    weirdness_measures$cross_index <- sapply(dfs, function(i) {
      # lapply(seq_len(ncol(i))[-1], function(z) {
      ss <- sapply(colnames(i)[-1], function(z) weirdness(i,
                                                          expected_sign = expected_sign,
                                                          response = z,
                                                          measure = "where_crosses_index"),
                   simplify = FALSE, USE.NAMES = TRUE)
      if(wmean) ss[[response]] else ss
    }, simplify = FALSE, USE.NAMES = TRUE)

    weirdness_measures$where_crosses <- sapply(zoi_vars_unique, function(i) {
      sapply(seq_along(weirdness_measures$cross_index[[i]]), function(z) {
        dfs[[i]][[1]][weirdness_measures$cross_index[[i]][[z]]]
      })#, simplify = FALSE, USE.NAMES = TRUE)
    }, simplify = FALSE, USE.NAMES = TRUE)

    weirdness_measures$n_crosses <- sapply(dfs, function(i) {
      # sapply(seq_along(weirdness_measures$cross_index[[i]]), function(z) {
      #   dfs[[i]][[1]][weirdness_measures$cross_index[[i]][[z]]]
      # })
      ss <- sapply(colnames(i)[-1], function(z) weirdness(i,
                                                          expected_sign = expected_sign,
                                                          response = z,
                                                          measure = "n_crosses"))
      if(wmean) ss[names(ss) == response][[1]] else ss
    })

    weirdness_measures$n_crosses_total <- sum(weirdness_measures$n_crosses)
    # }

  }

  # compute area
  if(any(grepl("response_area_opposite|response_area_ratio", measure))) {

    if(wmean) {
      weirdness_measures$response_area_opposite <- sapply(dfs, function(i) {
        weirdness(i,
                  expected_sign = expected_sign,
                  response = response,
                  measure = "response_area_opposite")
      })
      weirdness_measures$response_area_opposite_total <- sum(weirdness_measures$response_area_opposite)

      weirdness_measures$response_area_ratio <- sapply(dfs, function(i) {
        weirdness(i,
                  expected_sign = expected_sign,
                  response = response,
                  measure = "response_area_ratio")
      })
      weirdness_measures$response_area_ratio_total <- sum(weirdness_measures$response_area_ratio)
    } else {
      # implement later
      weirdness_measures$response_area_opposite <- paste0("This needs to be implemented for individual models. Please raise na issue on our Github repo.")
      weirdness_measures$response_area_opposite_total <- weirdness_measures$response_area_ratio <- weirdness_measures$response_area_ratio_total <- weirdness_measures$response_area_opposite
    }

  }

  if(any(grepl("n_inflection", measure))) {
    if(wmean) {
      weirdness_measures$n_inflection <- sapply(dfs, function(i) {
        weirdness(i,
                  expected_sign = expected_sign,
                  response = response,
                  measure = "n_inflection")
      })
      weirdness_measures$n_inflection_total <- sum(weirdness_measures$n_inflection)
    } else {
      weirdness_measures$n_inflection <- sapply(dfs, function(i) {
        sapply(colnames(i)[-1], function(z) weirdness(i,
                                                      expected_sign = expected_sign,
                                                      response = z,
                                                      measure = "n_crosses"))
      }) |>
        apply(MARGIN = 2, FUN = get(which_n_cross))
      weirdness_measures$n_inflection_total <- sum(weirdness_measures$n_inflection)
    }

    if(any(grepl("difference_inflection", measure))) {
      if(wmean) {
        weirdness_measures$difference_inflection <- sapply(dfs, function(i) {
          weirdness(i,
                    expected_sign = expected_sign,
                    response = response,
                    measure = "difference_inflection")
        })
        weirdness_measures$difference_inflection_total <- sum(weirdness_measures$difference_inflection)
      } else {
        # THIS NEEDS TO BE FURTHER IMPLEMENTED IF WE WANT TO IDENTIFY
        # THERE THE INFLECTION OCCURS
        weirdness_measures$difference_inflection <- sapply(dfs, function(i) {
          sapply(colnames(i)[-1], function(z) weirdness(i,
                                                        expected_sign = -1,
                                                        response = z,
                                                        measure = "difference_inflection"))
        }) |>
          apply(MARGIN = 2, FUN = get(which_n_cross))
        weirdness_measures$difference_inflection_total <- sum(weirdness_measures$difference_inflection)
      }
    }
  }

  # return
  weirdness_measures
}

#' Find inflection points in a curve
#'
#' Helper function to find inflection points in a curve, using the response axis values.
#'
#' @param `x` `[vector,numeric]` \cr Vector of `y` values, i.e., the response variable of
#' a curve.
#'
#' @keywords internal
#' @export
inflection <- function(x) c(FALSE, diff(diff(x) > 0) != 0)
