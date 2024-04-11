#' Predict bag of models in space
#'
#' @param gid `[string(1)="gid"]` \cr String with the name of the "gid" or point ID
#' column in `data`.
#'
#' @export
bag_predict_spat <- function(bag,
                             data,
                             # model = c("suit", "perm")[1],
                             input_type = c("df", "rast")[1], # only df implemented
                             output_type = c("df", "rast")[2],
                             baseline_max = NULL,
                             gridalign = TRUE,
                             gid = "gid",
                             coords = c("x33", "y33"),
                             crs = NULL,
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
  indiv_preds <- predvals * rep(bag$weights[good_models], each = nrow(predvals))
  # average (sum of individual weighted models)
  grd$linpred_ind_avg <- apply(indiv_preds, 1, sum, na.rm = TRUE)
  # uncertainty (sd of individual weighted models)
  grd$linpred_ind_sd <- apply(indiv_preds, 1, sd, na.rm = TRUE)
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
  out <- list(grid = grd, baseline_max = baseline_max,
              r_weighted_avg_pred = NULL,
              r_ind_summ_pred = NULL,
              r_ind_pred = NULL)

  #---
  # rasterize
  if(output_type == "rast"){

    if(verbose) print("Rasterizing...")

    # align pixels
    if (gridalign) { grd[, coords] <- grd[, coords]-50 }

    #---
    # average weighted model
    # rasterize
    grd_rast_wavg <- terra::rast(grd[, c(coords, "linpred_wavg")], type = "xyz", crs = crs)
    # exp values
    grd_rast_wavg <- exp(grd_rast_wavg)
    # quantile(exp(grd$linpred), prob = 0.975, na.rm = T)
    # truncate at 97.5% quantile of the baseline scenario
    if(is.null(baseline_max)) {
      baseline_max <- as.numeric(terra::global(grd_rast_wavg, fun = quantile, prob = 0.975, na.rm = T))
      # out$baseline_max <- baseline_max
    }
    grd_rast_wavg <- grd_rast_wavg/baseline_max
    grd_rast_wavg <- terra::ifel(grd_rast_wavg > 1, 1, grd_rast_wavg)
    grd_rast_wavg <- raster_rescale(grd_rast_wavg, to = c(0, 1))
    names(grd_rast_wavg) <- "suit_exp_weighted_avg"
    out$r_weighted_avg_pred <- grd_rast_wavg

    #---
    # summary of individual weighted models
    ind_summs <- c("ind_avg", "ind_sd")
    ind_vars <- paste0("linpred_", ind_summs)
    ind_names <- paste0("suit_exp_", ind_summs)
    i <- 1
    for(i in seq_along(ind_summs)) {

      # rasterize
      grd_rast <- terra::rast(grd[, c(coords, ind_vars[i])], type = "xyz", crs = crs)

      # exp values
      grd_rast <- exp(grd_rast)
      # quantile(exp(grd$linpred), prob = 0.975, na.rm = T)
      # truncate at 97.5% quantile of the baseline scenario
      if(is.null(baseline_max)) {
        baseline_max <- as.numeric(terra::global(grd_rast, fun = quantile, prob = 0.975, na.rm = T))
        # out$baseline_max <- baseline_max
      }
      grd_rast <- grd_rast/baseline_max
      grd_rast <- terra::ifel(grd_rast > 1, 1, grd_rast)
      grd_rast <- raster_rescale(grd_rast, to = c(0, 1))
      names(grd_rast) <- ind_names[i]

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
    i <- 1
    for(i in seq_along(ind_summs)) {

      # rasterize
      grd_rast <- terra::rast(grd[, c(coords, ind_vars[i])], type = "xyz", crs = crs)

      # exp values
      grd_rast <- exp(grd_rast)
      # quantile(exp(grd$linpred), prob = 0.975, na.rm = T)
      # truncate at 97.5% quantile of the baseline scenario
      if(is.null(baseline_max)) {
        baseline_max <- as.numeric(terra::global(grd_rast, fun = quantile, prob = 0.975, na.rm = T))
        # out$baseline_max <- baseline_max
      }
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
