#' @name influence_functions
#'
#' @title Influence functions
#'
#' @description Decay influence functions, parameterized based on the zone of influence (ZoI).
#'
#' @example examples/influence_functions_example.R
#'
#' @rdname influence_functions
#' @export
exp_decay <- function(x, zoi = NULL, exp_decay_parms = c(1, 0.01),
                      zoi_decay_threshold = 0.05,
                      origin = 0,
                      half_life = NULL, zoi_hl_ratio = NULL,
                      oneside = TRUE) {

  # define lambda depending on the input parameter
  if(!is.null(zoi)) {

    if(is.null(zoi_hl_ratio)) {
      lambda <- log(1/zoi_decay_threshold) / zoi
    } else {
      half_life <- zoi/zoi_hl_ratio
      lambda <- log(2)/half_life
    }

  } else {

    if(!is.null(half_life)) {
      lambda <- log(2)/half_life
    } else {
      lambda <- exp_decay_parms[2]
    }
  }

  # return function
  if(oneside) func <- identity else func <- abs

  exp_decay_parms[1] * exp(- lambda * func(x - origin))
}

#' @rdname influence_functions
#' @export
threshold_decay <- function(x, zoi, constant_influence = 1, origin = 0, oneside = TRUE) {
  if(oneside) func <- identity else func <- abs
  ifelse(func(x - origin) < zoi, constant_influence, 0)
}

#' @rdname influence_functions
#' @export
step_decay <- threshold_decay

#' @rdname influence_functions
#' @export
bartlett_decay <- function(x, zoi, intercept = 1, origin = 0, oneside = TRUE) {
  if(oneside) func <- identity else func <- abs
  beta = -intercept/zoi
  ifelse(func(x - origin) < zoi, intercept + beta * func(x - origin), 0)
}

#' @rdname influence_functions
#' @export
tent_decay <- bartlett_decay

#' @rdname influence_functions
#' @export
linear_decay <- bartlett_decay

#' @rdname influence_functions
#' @export
gaussian_decay <- function(x, zoi = NULL, hnorm_decay_parms = c(1, 0.01),
                           sigma = NULL, zoi_decay_threshold = 0.05,
                           origin = 0, ...) {

  if(!is.null(zoi)) {
    lambda = log(1/zoi_decay_threshold) / (zoi**2)
  } else {
    if(!is.null(sigma)) {
      lambda = 1/(2*sigma**2)
    } else {
      lambda <- hnorm_decay_parms[2]
    }

  }

  hnorm_decay_parms[1] * exp(- lambda * (x - origin)**2)
}

#' @rdname influence_functions
#' @export
half_norm_decay <- gaussian_decay
