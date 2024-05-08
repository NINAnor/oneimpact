#' Predict bag of models in space
#'
#' @param bag `[bag,list]` \cr A bag of models, resulting from a call to [oneimpact::bag_models()].
#' @param data `[data.frame,SpatRaster]` \cr The spatial grid or stack of rasters
#' with the layers in space, to be used for the spatial prediction. All variables in the
#' bag$formula must be present in `data`.
#' @param input_type `[string(1)="df"]{"df","rast"}` \cr Type of input object. Either a
#' `data.frame` with the predictor variables as columns (if `input_type = "df"`, default),
#' or a `SpatRaster` with the predictor variables as layers (if `input_type = "rast"`).
#' So far, only `input_type = "df"` is implemented.
#' @param output_type `[string(1)="rast"]{"df","rast"}` \cr Type of output object.
#' Typically, the same type of object as `input_type`, but rasters can also be saved
#' if the input is a `data.frame` when `output_type = "rast"`, and data.frames can
#' also be saved if the input is a `SpatRaster` when `output_type = "df"`.
#' @param gid `[string(1)="gid"]` \cr String with the name of the "gid" or point ID
#' column in `data`. Only relevant if `input_type = "df"`.
#' @param coords `[vector,string(2)=c("x", "y")]` \cr Vector with two elements with
#' the names of the coordinates representing (x,y) coordinates for the pixels.
#' Only relevant if `input_type = "df"`.
#' @param crs `[string(1)=NULL]` \cr Code for the coordinate reference system of the
#' output raster. Only relevant if `input_type = "df"`. For more details, check
#' [terra::crs()].
#'
#'
#' @export
bag_predict_spat <- function(bag,
                             data,
                             # model = c("suit", "perm")[1],
                             input_type = c("df", "rast")[1], # only df implemented
                             output_type = c("df", "rast")[2],
                             gridalign = TRUE,
                             gid = "gid",
                             coords = c("x33", "y33"),
                             crs = NULL,
                             prediction_max_quantile = 0.999,
                             uncertainty_quantiles = c(0.25, 0.75),
                             plotit = FALSE,
                             verbose = FALSE) {

  #---
  # prepare

  if(verbose) print("Preparing data...")

  # model covariates and response
  model_vars <- all.vars(bag$formula_no_strata)[-1]
  case <- extract_response_strata(bag$formula)$response

  # get data for prediction
  # cols %in% names(data)
  grd <- data[, c(gid, coords, cols)]
  grd[[case]] <- rep(0, nrow(grd))
  nrow(grd)
  colnames(grd)

  nrow(grd <- na.omit(grd))

  #---
  # prepare model matrix
  mm <- model.matrix(bag$formula_no_strata, data = grd)
  nrow(mm)

  #---
  # predict

  if(verbose) print("Predicting...")

  # predict only for models with weight > 0
  good_models <- which(bag$weights > 0)
  # coefs
  coefs <- bag$coef[, good_models, drop = FALSE]
  # prediction
  predvals <- (mm %*% coefs)

  # average weighted model prediction
  grd$linpred_wavg <- predvals %*% bag$weights[good_models]

  # prediction for each model
  # checking it is right
  # all(rep(bag_summary$weights[good_models], each = nrow(predvals)) == as.vector(matrix(rep(bag_summary$weights[good_models], each = nrow(predvals)), nrow = nrow(predvals))))
  # indiv_preds <- predvals * rep(bag$weights[good_models], each = nrow(predvals))
  indiv_preds <- predvals
  # average (sum of individual weighted models)
  # grd$linpred_ind_w_med <- apply(indiv_preds, 1, modi::weighted.quantile, w = bag$weights[good_models],
  #                                probs = 0.5)
  grd$linpred_ind_w_med <- apply(indiv_preds, 1, DescTools::Quantile, weights = bag$weights[good_models],
                                 type = 5, probs = 0.5)
  # uncertainty (sd of individual weighted models)
  # grd$linpred_ind_w_quart_min <- apply(indiv_preds, 1, modi::weighted.quantile, w = bag$weights[good_models],
  #                                    prob = uncertainty_quantiles[1])
  grd$linpred_ind_w_quart_min <- apply(indiv_preds, 1, DescTools::Quantile, weights = bag$weights[good_models],
                                 type = 5, probs = uncertainty_quantiles[1])
  # grd$linpred_ind_w_quart_max <- apply(indiv_preds, 1, modi::weighted.quantile, w = bag$weights[good_models],
  #                                    prob = uncertainty_quantiles[2])
  grd$linpred_ind_w_quart_max <- apply(indiv_preds, 1, DescTools::Quantile, weights = bag$weights[good_models],
                                       type = 5, probs = uncertainty_quantiles[2])

  # individual predictions
  new_cols <- ncol(grd)+1:ncol(indiv_preds)
  grd[new_cols] <- indiv_preds
  indiv_names <- paste0("linpred_ind_", names(good_models))
  names(grd)[new_cols] <- indiv_names
  names(grd)

  #---
  # come back to space

  #---
  # output object
  out <- list(grid = grd,
              weights = bag$weights[good_models],
              r_weighted_avg_pred = NULL,
              r_ind_summ_pred = NULL,
              r_ind_pred = NULL)

  #---
  # rasterize
  if(output_type == "rast"){

    if(verbose) print("Rasterizing...")

    # align pixels
    # if (gridalign) { grd[, coords] <- grd[, coords]-50 }

    #---
    # average weighted model
    # rasterize
    grd_rast_wavg <- terra::rast(grd[, c(coords, "linpred_wavg")], type = "xyz", crs = crs)
    # exp values
    grd_rast_wavg <- exp(grd_rast_wavg)
    # quantile(exp(grd$linpred), prob = 0.975, na.rm = T)
    # truncate at 97.5% quantile of the baseline scenario
    baseline_max <- as.numeric(terra::global(grd_rast_wavg, fun = quantile, prob = prediction_max_quantile, na.rm = T))
    grd_rast_wavg <- grd_rast_wavg/baseline_max
    grd_rast_wavg <- terra::ifel(grd_rast_wavg > 1, 1, grd_rast_wavg)
    grd_rast_wavg <- raster_rescale(grd_rast_wavg, to = c(0, 1))
    names(grd_rast_wavg) <- "suit_exp_weighted_avg"
    out$r_weighted_avg_pred <- grd_rast_wavg
    # plot(grd_rast_wavg)

    #---
    # summary of individual weighted models
    ind_summs <- c("ind_w_med", "ind_w_iqr")
    ind_vars <- paste0("linpred_", c("ind_w_med", c("ind_w_quart_min", "ind_w_quart_max")))
    ind_names <- paste0("suit_", c("exp_", ""), ind_summs)
    # i <- 2
    for(i in seq_along(ind_summs)) {

      # rasterize
      if(i == 1) {
        grd_rast <- terra::rast(grd[, c(coords, ind_vars[i])], type = "xyz", crs = crs)
        # exp values
        grd_rast <- exp(grd_rast)
      } else {
        grd_rast <- terra::rast(grd[, c(coords, ind_vars[c(2,3)])], type = "xyz", crs = crs)
        # diff
        grd_rast <- terra::diff(grd_rast)
        # exp values
        # grd_rast <- exp(grd_rast)
      }
      # plot(grd_rast)

      # quantile(exp(grd$linpred), prob = 0.975, na.rm = T)
      # truncate at 97.5% quantile of the baseline scenario
      baseline_max <- as.numeric(terra::global(grd_rast, fun = quantile, prob = prediction_max_quantile, na.rm = T))

      grd_rast <- grd_rast/baseline_max
      grd_rast <- terra::ifel(grd_rast > 1, 1, grd_rast)
      grd_rast <- raster_rescale(grd_rast, to = c(0, 1))
      names(grd_rast) <- ind_names[i]
      # plot(grd_rast)

      if(i == 1) {
        # rasterize
        out$r_ind_summ_pred <- grd_rast
      } else {
        # rasterize and concatenate
        out$r_ind_summ_pred <- c(out$r_ind_summ_pred, grd_rast)
      }
    }

    #---
    # individual predictions
    ind_vars <- grep("linpred_ind_Resample", names(grd), value = TRUE)
    ind_summs <- sub("linpred_", "", ind_vars)
    ind_names <- paste0("suit_exp_", ind_summs)
    # i <- 1
    for(i in seq_along(ind_summs)) {

      # rasterize
      grd_rast <- terra::rast(grd[, c(coords, ind_vars[i])], type = "xyz", crs = crs)

      # exp values
      grd_rast <- exp(grd_rast)
      # quantile(exp(grd$linpred), prob = 0.975, na.rm = T)
      # truncate at 97.5% quantile of the baseline scenario
      baseline_max <- as.numeric(terra::global(grd_rast, fun = quantile, prob = prediction_max_quantile, na.rm = T))

      grd_rast <- grd_rast/baseline_max
      grd_rast <- terra::ifel(grd_rast > 1, 1, grd_rast)
      grd_rast <- raster_rescale(grd_rast, to = c(0, 1))
      names(grd_rast) <- ind_names[i]

      if(i == 1) {
        # rasterize
        out$r_ind_pred <- grd_rast
      } else {
        # rasterize and concatenate
        out$r_ind_pred <- c(out$r_ind_pred, grd_rast)
      }

    }

  }

  if(plotit) plot(out$r_weighted_avg_pred)
  if(plotit) plot(out$r_ind_summ_pred)
  if(plotit) plot(out$r_ind_pred)

  # return
  out
}
