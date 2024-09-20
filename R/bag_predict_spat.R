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
                             prediction_type = c("exp", "exponential", "linear")[1],
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

  if(input_type == "rast") {

    # get crs
    if(is.null(crs)) {
      crs <- terra::crs(data)
    } else {
      if ((crs == terra::crs('epsg:25833')) == FALSE) {
        warning("The argument 'crs' is not the same as the CRS for you 'data' raster. ")
      }
    }

    data <- terra::as.data.frame(data, cells = TRUE, xy = TRUE, na.rm = FALSE)
    names(data)[1:3] <- c(gid, coords[1], coords[2])


  }

  # get data for prediction
  # cols %in% names(data)
  grd <- data[, c(gid, coords, model_vars)]
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
  if(grepl("exp", prediction_type)) grd$linpred_wavg <- exp(grd$linpred_wavg)

  # prediction for each model
  # checking it is right
  # all(rep(bag_summary$weights[good_models], each = nrow(predvals)) == as.vector(matrix(rep(bag_summary$weights[good_models], each = nrow(predvals)), nrow = nrow(predvals))))
  # indiv_preds <- predvals * rep(bag$weights[good_models], each = nrow(predvals))
  indiv_preds <- predvals
  if(grepl("exp", prediction_type)) indiv_preds <- exp(indiv_preds) ## NEW
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
    # individual predictions
    ind_vars <- grep("linpred_ind_Resample", names(grd), value = TRUE)
    ind_summs <- sub("linpred_", "", ind_vars)
    ind_names <- paste0("suit_exp_", ind_summs)
    # i <- 1
    for(i in seq_along(ind_summs)) {

      # rasterize
      grd_rast <- terra::rast(grd[, c(coords, ind_vars[i])], type = "xyz", crs = crs)

      # exp values
      # grd_rast <- exp(grd_rast)
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

    #---
    # summary of individual weighted models
    ind_summs <- c("ind_w_med", "ind_w_iqr", "ind_w_iqr_old")
    ind_vars <- paste0("linpred_", c("ind_w_med", c("ind_w_quart_min", "ind_w_quart_max")))
    ind_names <- paste0("suit_", c("exp_", ""), ind_summs)
    # i <- 2
    for(i in seq_along(ind_summs)) {

      # rasterize
      if(i == 1) {
        grd_rast <- terra::rast(grd[, c(coords, ind_vars[i])], type = "xyz", crs = crs)
      } else {
        if(i == 3) {
          grd_rast <- terra::rast(grd[, c(coords, ind_vars[c(2,3)])], type = "xyz", crs = crs)
        } else {
          grd_rast <- app(out$r_ind_pred, function(x)
            DescTools::Quantile(x, weights = bag$weights[good_models],
                                type = 5, probs = uncertainty_quantiles))
        }
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


  }

  # plot(out$r_ind_summ_pred[[2]])
  # plot(terra::diff(app(out$r_ind_pred, function(x) terra::quantile(x, probs = c(0.05, 0.95), na.rm = TRUE))))
  # plot(app(out$r_ind_pred, sd))
  # plot(out$r_ind_pred[[grepl("068|086|047|020|133|100", names(out$r_ind_pred))]])

  if(plotit) plot(out$r_weighted_avg_pred)
  if(plotit) plot(out$r_ind_summ_pred)
  if(plotit) plot(out$r_ind_pred)

  # return
  out
}

#' @export
bag_predict_spat_vars <- function(bag,
                                  data,
                                  predictor_table_zoi,
                                  # model = c("suit", "perm")[1],
                                  input_type = c("df", "rast")[1], # only df implemented
                                  output_type = c("df", "rast")[2],
                                  prediction_type = c("exp", "exponential", "linear")[1],
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
  grd <- data[, c(gid, coords, model_vars)]
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

  #---
  # output object
  out <- list(vars = NULL,
              grid = NULL,
              weights = bag$weights[good_models],
              r_weighted_avg_pred = NULL,
              r_ind_summ_pred = NULL,
              r_ind_pred = NULL)

  # lop for variables, make all others zero
  vars <- unique(predictor_table_zoi$variable)
  for(i in seq_along(vars)) {

    if(verbose) print(paste0("Computing maps for variable ", vars[i], "..."))

    # get variables
    var_i <- vars[i]
    vars_full <- predictor_table_zoi$term_zoi[predictor_table_zoi$variable == var_i]

    # copy table
    coefs1 <- coefs
    grd1 <- grd

    if(length(vars_full) > 0) {

      rows_i <- which(grepl(paste0(vars_full, collapse = "|"), rownames(coefs1)))
      # other rows
      rows_non_i <- which(!grepl(paste0(vars_full, collapse = "|"), rownames(coefs1)))
      # make them zero
      coefs1[rows_non_i,] <- 0

      # prediction
      predvals <- (mm %*% coefs1)

      # average weighted model prediction
      col_name <- paste0("pred_", var_i, "wavg")
      grd1[col_name] <- predvals %*% bag$weights[good_models]
      if(grepl("exp", prediction_type)) grd1[col_name] <- exp(grd1[col_name])

      # prediction for each model
      indiv_preds <- predvals
      if(grepl("exp", prediction_type)) indiv_preds <- exp(indiv_preds)

      # weighted median
      col_name <- paste0("pred_", var_i, "w_med")
      grd1[col_name] <- apply(indiv_preds, 1, DescTools::Quantile, weights = bag$weights[good_models],
                              type = 5, probs = 0.5)

      # individual predictions
      new_cols <- ncol(grd1)+1:ncol(indiv_preds)
      grd1[new_cols] <- indiv_preds
      indiv_names <- paste0("pred_", var_i, "ind_", names(good_models))
      names(grd1)[new_cols] <- indiv_names
      names(grd1)

      #---
      # come back to space

      out$vars[[i]] <- var_i
      out$grid[[i]] <- grd1[, c(gid, coords)]
      out$grid[[i]] <- cbind(out$grid[[i]],
                             grd1[, grep(paste0("pred_", var_i), colnames(grd1))])

      #---
      # rasterize
      if(output_type == "rast"){

        if(verbose) print("Rasterizing...")

        # align pixels
        # if (gridalign) { grd[, coords] <- grd[, coords]-50 }

        #---
        # average weighted model
        # rasterize
        col_name <- paste0("pred_", var_i, "wavg")
        grd_rast_wavg <- terra::rast(grd1[, c(coords, col_name)], type = "xyz", crs = crs)
        # quantile(exp(grd$linpred), prob = 0.975, na.rm = T)
        # truncate at 97.5% quantile of the baseline scenario
        # baseline_max <- as.numeric(terra::global(grd_rast_wavg, fun = quantile, prob = prediction_max_quantile, na.rm = T))
        # grd_rast_wavg <- grd_rast_wavg/baseline_max
        # grd_rast_wavg <- terra::ifel(grd_rast_wavg > 1, 1, grd_rast_wavg)
        # grd_rast_wavg <- raster_rescale(grd_rast_wavg, to = c(0, 1))
        names(grd_rast_wavg) <- paste0("suit_exp_", var_i, "wavg")
        out$r_weighted_avg_pred[[i]] <- grd_rast_wavg
        # plot(grd_rast_wavg)

        # weighted median
        # col_name <- paste0("pred_", var_i, "w_med")
        # grd_rast_w_med <- terra::rast(grd1[, c(coords, col_name)], type = "xyz", crs = crs)
        # names(grd_rast_w_med) <- paste0("suit_exp_", var_i, "wmed")
        # out$r_ind_summ_pred[[i]] <- grd_rast_w_med
        # # plot(grd_rast_w_med)

        #---
        # individual predictions
        ind_vars <- grep(paste0("pred_", var_i, "ind_Resample"), names(grd1), value = TRUE)
        ind_summs <- sub("pred_", "", ind_vars)
        ind_names <- paste0("suit_exp_", ind_summs)
        # i <- 1
        for(j in seq_along(ind_summs)) {

          # rasterize
          grd_rast <- terra::rast(grd1[, c(coords, ind_vars[j])], type = "xyz", crs = crs)

          # exp values
          # grd_rast <- exp(grd_rast)
          # quantile(exp(grd$linpred), prob = 0.975, na.rm = T)
          # truncate at 97.5% quantile of the baseline scenario
          # baseline_max <- as.numeric(terra::global(grd_rast, fun = quantile, prob = prediction_max_quantile, na.rm = T))

          # grd_rast <- grd_rast/baseline_max
          # grd_rast <- terra::ifel(grd_rast > 1, 1, grd_rast)
          # grd_rast <- raster_rescale(grd_rast, to = c(0, 1))
          names(grd_rast) <- ind_names[j]
          # plot(grd_rast)

          if(j == 1) {
            # rasterize
            out$r_ind_pred[[i]] <- grd_rast
          } else {
            # rasterize and concatenate
            out$r_ind_pred[[i]] <- c(out$r_ind_pred[[i]], grd_rast)
          }

        }

        # #---
        # summary of individual weighted models
        ind_summs <- paste0("pred_", var_i, c("w_med", "ind_w_iqr"))
        ind_vars <- ind_summs
        ind_names <- paste0("suit_exp_", sub("pred_", "", ind_vars))
        # i <- 2
        for(j in seq_along(ind_summs)) {

          # rasterize
          if(j == 1) {
            grd_rast <- terra::rast(grd1[, c(coords, ind_vars[j])], type = "xyz", crs = crs)
          } else {
            grd_rast <- app(out$r_ind_pred[[i]], function(x)
              DescTools::Quantile(x, weights = bag$weights[good_models],
                                  type = 5, probs = uncertainty_quantiles))

            # diff
            grd_rast <- terra::diff(grd_rast)
            # exp values
            # grd_rast <- exp(grd_rast)
          }
          # plot(grd_rast)

          # quantile(exp(grd$linpred), prob = 0.975, na.rm = T)
          # truncate at 97.5% quantile of the baseline scenario
          # baseline_max <- as.numeric(terra::global(grd_rast, fun = quantile, prob = prediction_max_quantile, na.rm = T))
          # grd_rast <- grd_rast/baseline_max
          # grd_rast <- terra::ifel(grd_rast > 1, 1, grd_rast)
          # grd_rast <- raster_rescale(grd_rast, to = c(0, 1))
          names(grd_rast) <- ind_names[j]
          # plot(grd_rast)

          if(j == 1) {
            # rasterize
            out$r_ind_summ_pred[[i]] <- grd_rast
          } else {
            # rasterize and concatenate
            out$r_ind_summ_pred[[i]] <- c(out$r_ind_summ_pred[[i]], grd_rast)
          }
        }


      }

    }


  }

  out

}
