#' Get estimates of ZOI from response curves
#'
#' @example examples/zoi_from_curve_example.R
#'
#' @export
zoi_from_curve <- function(x,
                           percentage = 0.95,
                           curve = c("median", "mean"),
                           ci = TRUE,
                           type = c("linear", "exp")[1],
                           # n_features = 1,
                           mean_col_name = "mean",
                           median_col_name = "quantile:0.5",
                           ci_col_name = c("quantile:0.025", "quantile:0.975")) {
  UseMethod("zoi_from_curve")
}

#' @export
zoi_from_curve <- function(x,
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

    # max effect size
    max_effect_size_main <- main_response[which.max(abs(main_response))]

    # zoi radius main
    y_value_percentage <- (1 - percentage) * (max_effect_size_main)
    x_radius_index <- which(abs(main_response) < abs(y_value_percentage))[1] - 1
    zoi_radius_main <- x[[xvar]][x_radius_index]

    # impact
    x_vals <- x[[xvar]][1:x_radius_index]
    y_vals <- main_response[1:x_radius_index]
    if(type == "exp") y_vals <- exp(y_vals) - 1
    y_vals <- abs(y_vals) - min(abs(y_vals)) # get positive and discount are above y(ZOI)

    impact_main <- DescTools::AUC(x_vals, y_vals)

    max_effect_size[id] <- ifelse(type == "exp", exp(max_effect_size_main), max_effect_size_main)
    zoi_radius[id] <- zoi_radius_main
    effect_zoi_radius[id] <- ifelse(type == "exp", exp(y_value_percentage), y_value_percentage)
    impact[id] <- impact_main
  }

  # repeat that for the CI
  if(ci) {

    ci_id <- 1

    for(ci_id in seq(ci_col_name)) {

      ci_response <- x[[ci_col_name[ci_id]]]

      # max effect size
      ci_max <- ci_response[which.max(abs(main_response))]

      # zoi radius main
      # y_value_percentage <- (1 - percentage) * max_effect_size
      x_radius_ci_index <- which(abs(ci_response) < abs(y_value_percentage))[1] - 1
      zoi_radius_ci <- x[[xvar]][x_radius_ci_index]

      # impact
      x_vals <- x[[xvar]][1:x_radius_ci_index]
      y_vals <- ci_response[1:x_radius_ci_index]
      if(type == "exp") y_vals <- exp(y_vals) - 1
      y_vals <- abs(y_vals) - min(abs(y_vals)) # get positive and discount are above y(ZOI)

      impact_ci <- DescTools::AUC(x_vals, y_vals)

      max_effect_size[ci_id + 2] <- ifelse(type == "exp", exp(ci_max), ci_max)
      zoi_radius[ci_id + 2] <- zoi_radius_ci
      effect_zoi_radius[ci_id + 2] <- ifelse(type == "exp", exp(y_value_percentage), y_value_percentage)
      impact[ci_id + 2] <- impact_ci
    }
  }

  out <- list(max_effect_size = max_effect_size,
              zoi_radius = zoi_radius,
              effect_zoi_radius = effect_zoi_radius,
              impact = impact) |>
    dplyr::bind_rows() |>
    t()

  colnames(out) <- c("mean", "median", ci_col_name[1], ci_col_name[2])
  out
}


# test n_features
# test exp vs linear
# implement function.df this one
# implement function var
# implement function bag - all vars


