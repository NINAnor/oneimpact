#' @name zoi_functions
#'
#' @title Zone of Influence (ZOI) functions
#'
#' @description Computes the decay functions that represent the Zone of
#' Influence (ZOI). The functions' radius (parameter `radius`)
#' controls how far the zone of influence of an infrastructure or disturbance reaches,
#' and the functions' shape (parameter `type`) represent represent how the ZOI
#' decays in space. Given a function shape (`type`) is chosen, the rate
#' of decay of the different ZOI functions is parameterized based on
#' the ZOI radius -- e.g the slope of [oneimpact::linear_decay()] is defined
#' so that the function decreases to zero at the ZOI radius.
#' These functions can be used to transform arrays of (Euclidean)
#' distance values (in one dimension) or rasters of (Euclidean) distance
#' (in two dimensions) into ZOI values. The distances might represent
#' the distance to human infrastructure, sources of disturbance, or
#' more broadly any type of land use class or spatial variable.
#'
#' @details
#' A generic function [oneimpact::dist_decay()] can be used to compute
#' ZOI values according to functions with different shapes (parameter `type`)
#' and radii (parameter `radius`). Alternatively, there are specific functions
#' implemented for each ZOI shape.
#'
#' For the threshold function ([oneimpact::threshold_decay()]) and the linear decay
#' function ([oneimpact::linear_decay()]), the ZOI radius (parameter `radius`) is the
#' distance `x` where the ZOI function value decreases to zero.
#' For the linear decay, this is done by setting
#' the slope of the linear function as `-intercept/radius`, where `intercept`
#' is the intercept of the linear function (here, the maximum value at `x = 0`).
#'
#' For non-vanishing functions that approach zero asymptotically
#' ([oneimpact::exp_decay()], [oneimpact::gaussian_decay()]), a certain limit value must be given to define
#' the ZOI radius -- so that the ZOI radius is defined as the distance `x` where the
#' ZOI function goes below this limit value. For these functions,
#' different parameters are available
#' for setting the relationship between the ZOI function value and the ZOI radius.
#'
#' Some functions have multiple possible names, for the sake of flexibility:
#' - [oneimpact::linear_decay()], [oneimpact::bartlett_decay()], and
#' [oneimpact::tent_decay()] are the same function;
#' - [oneimpact::threshold_decay()] and [oneimpact::step_decay()] are the same function;
#' - [oneimpact::gaussian_decay()] and [oneimpact::half_norm_decay()] are the same function.
#'
#' Alternatively, [oneimpact::dist_decay()] can call all of them, given a
#' ZOI shape is specified through the parameter `type`.
#'
#' Other functions might be implemented in the future.
#'
#' # Definitions
#'
#' Here are some formal definitions for the ZOI functions \eqn{\phi(d_i, r)},
#' where \eqn{d_i} is the distance to the feature \eqn{i} of an infrastructure or
#' source of disturbance and \eqn{r} is the radius of the zone of influence:
#' - `threshold_decay()`: the threshold or step decay function \eqn{\phi_{threshold}} is
#' positive and constant within the ZOI radius \eqn{r}, and null for \eqn{x \ge r}:
#' \deqn{
#' \phi_{threshold}(d_i, r_k) = c \text{ if } d_i < r, 0 \text{ otherwise}
#' }
#' where \eqn{c} is a constant value (by default `c = 1`).
#' - `linear_decay()`: the linear (or tent/Bartlett) decay function \eqn{\phi_{linear}}
#' decreases
#' linearly from a maximum value \eqn{c} (the intercept, by default `c = 1`) to
#' zero when \eqn{x \ge r}:
#' \deqn{\phi_{linear}(d_i, r) = c - c/r \text{ if } x < r, \text{ 0 otherwise}}
#' - `exp_decay()`: the exponential decay function \eqn{\phi_{exp}} decreases
#' exponentially from a maximum value \eqn{c} (by default `c = 1`) with a rate
#' \eqn{\lambda}, which is defined by \eqn{r} and a ZOI limit value
#' \eqn{\phi_{lim}}, a small ZOI value below which the influence is considered negligible:
#' \deqn{\phi_{exp}(d_i, r, \phi_{lim}) = c exp(-\lambda d_i)}
#' with
#' \deqn{\lambda = ln(1/\phi_{lim}) / r}
#' In this context, the ZOI radius \eqn{r} is the distance beyond which
#' \eqn{\phi_{exp} < \phi_{lim}}.
#' - `gaussian_decay()`: the Gaussian decay function \eqn{\phi_{Gauss}}
#' follows a Gaussian (half-normal) decay with maximum \eqn{c} (by default `c = 1`)
#' and a decay rate \eqn{\lambda}
#' defined by \eqn{r} and a ZOI limit value
#' \eqn{\phi_{lim}}, a small ZOI value below which the influence is considered negligible:
#' \deqn{\phi_{Gauss}(d_i, r, \phi_{lim}) = c exp(-\lambda d_i^2)}
#' with
#' \deqn{\lambda = ln(1/\phi_{lim}) / (r^2)}
#' In this context, the ZOI radius \eqn{r} is the distance beyond which
#' \eqn{\phi_{exp} < \phi_{lim}}. Note that \eqn{\lambda} is defined differently
#' for the `gaussian_decay` and the `exp_decay` functions.
#'
#' # Parameterization
#'
#' Some of the shapes of the ZOI (parameter `type` in `dist_decay()`) might be
#' parameterized in multiple ways. Here is a brief description of each possibility:
#'
#' \itemize{
#'
#'   \item For the `"Gauss"` or `"half_norm"` shapes, the ZOI follows a half-normal shape: \cr
#'   `intercept * exp(-lambda * (euclidean_distance^2))`. `intercept` and `lambda` are
#'   parameters to be defined. There are three ways of specifying `lambda`:
#'   \itemize{
#'     \item If the `radius = NULL` (default), `lambda` is a parameter by itself to
#'     be specified by the user. In all other cases (below)
#'     the value of this parameter is ignored, even if provided.
#'     \item If the parameter `radius` is provided, the rate of decay is given by
#'     `lambda = log(1/zoi_limit) / (radius**2)`. In other words, `lambda` is defined
#'     so that the function decreases to `zoi_limit` when `x = radius`.
#'     \item If the `radius = NULL` and `sigma` is provided, `lambda` is defined as
#'     `lambda = 1/(2*sigma**2)`.
#'   }
#'
#'   \item For the `"exp_decay"` shape, the ZOI follows an exponential decay shape: \cr
#'   `intercept * exp(-lambda * euclidean_distance)`. `intercept` and `lambda` are
#'   parameters to be defined. There are four ways of specifying `lambda`:
#'   \itemize{
#'     \item If the `radius = NULL` (default), `lambda` is a parameter by itself to
#'     be specified by the user. In all other cases (below)
#'     the value of this parameter is ignored, even if provided.
#'     \item If the parameter `radius` is provided, the rate of decay is given by
#'     `lambda = log(1/zoi_limit) / radius`. In other words, `lambda` is defined
#'     so that the function decreases to `zoi_limit` when `x = radius`.
#'     \item If the `radius = NULL` and `half_life` is given, `lambda` is defined
#'     based on the half life of the exponential function -- the distance at which
#'     the function decreases to 1/2. If `zoi_hl_ratio = NULL`, `lambda` is defined
#'     as `lambda = log(2)/half_life`.
#'     \item The last possibility is to specify `zoi_hl_ratio`, the ratio between
#'     the ZOI radius and the half life of the exponential function. For instance,
#'     if `zoi_hl_ratio = 4`, this means the ZOI radius is defined as `4*half_life`.
#'     If `zoi_hl_ratio` is provided, the exponential `half_life` is defined based on
#'     this parameter and `lambda` is defined accordingly, based on the relationship
#'     above. In this case, `radius` is ignored, even if specified.
#'   }
#'
#'   \item For the `"bartlett"`, `"linear_decay"`, or `"tent_decay"` shapes, the ZOI follows a
#'   linear decay shape (`y = a*x + b`) within the ZOI radius (parameter `radius`).
#'   The intercept of the linear function (`b`) is given by the parameter `intercept`
#'   and the slope (`a`) is given by `-intercept/radius`.
#'
#'   \item For the `"threshold"` or `"step"` shapes, a constant influence is consider
#'   within the
#'   zone of influence radius (parameter `radius`). All pixels closer than
#'   `radius` to infrastructure are considered as "under the influence" of
#'   the nearest feature, with a constant influence value defined by the
#'   `intercept` parameter, and all values/pixels beyond `radius` are assumed to have
#'   zero influence.
#' }
#'
#' @param x `[numeric,SpatRaster,RasterLayer]` \cr Euclidean distance from an infrastructure, source
#' of disturbance, or feature/class of interest. It can be a single value, an array
#' of values, or a raster object. It must not necessarily be an Euclidean distance,
#' but preferably it should be a distance measured in meters, to ease interpretation
#' (e.g. geodesic distance).
#'
#' @param radius `[numeric(1)]` \cr Radius of the zone of influence (ZOI),
#' the distance at which the ZOI vanishes or goes below a given minimum limit value
#' `zoi_limit`. See details.
#'
#' @param zoi_limit `[numeric(1)=0.05]` \cr For non-vanishing functions
#' (e.g. `exp_decay`, `gaussian_decay`), this value is used to set the relationship
#' between the ZOI radius and the decay functions:
#' `radius` is defined as the minimum distance `x` at which the ZOI assumes values
#' below `zoi_limit`. The default is 0.05. This parameter is used only
#' if `radius` is not `NULL`.
#'
#' @param type `[character(1)="Gauss"]{"Gauss", "exp_decay", "bartlett",
#' "linear", "tent", "threshold", "step"}` \cr Type or shape of the decay distance.
#' \itemize{
#'   \item If `"Gauss"` or `"half_norm"`, the ZOI follows a half-normal shape: \cr
#'   `intercept * exp(-lambda * (euclidean_distance^2))`. `intercept` and `lambda` are
#'   parameters to be defined -- see details.
#'   \item If `"exp_decay"`, the ZOI follows an exponential decay shape: \cr
#'   `intercept * exp(-lambda * euclidean_distance)`. `intercept` and `lambda` are
#'   parameters to be defined -- see details.
#'   \item If `"bartlett"`, `"linear_decay"`, or `"tent_decay"`, the ZOI follows a
#'   linear decay shape within the ZOI radius (parameter `radius`).
#'   \item If `"threshold"` or `"step"`, a constant influence is consider within the
#'   zone of influence radius (parameter `radius`). All pixels closer than
#'   `radius` to infrastructure are considered as "under the influence" of
#'   the nearest feature, with a constant influence value defined by the
#'   `intercept` parameter, and all values/pixels beyond `radius` are
#'   assumed to have zero influence.
#' }
#'
#' @param intercept `[numeric(1)=1]` \cr Maximum value of the ZOI function at
#' when the distance from disturbance sources is zero (`x = 0`).
#' For the `threshold_decay` and `step_decay` functions, `intercept` is
#' the constant value of the Zone of Influence within the ZOI `radius`.
#' For the other ZOI functions, `intercept`
#' is the value of the functions at the origin (where the sources of disturbance
#' are located, i.e. `x = 0`).
#' Default is `intercept = 1`.
#'
#' @param origin `[numeric(1)=0]` \cr In which position (in 1 dimension) is located
#' the infrastructure or source of disturbance? Default is zero. For raster objects,
#' this parameter should be ignored.
#'
#' @param oneside `[logical(1)=TRUE]` \cr If `FALSE`, negative distance values
#' are considered symmetrically and their transformation is always positive.
#' This parameter is only meaningful if `x` is a vector of values, not a
#' raster object.
#'
#' @return The ZOI values for a given array of distance values if `x` is numeric,
#' or a raster object delimiting the ZOI if `x` corresponds to the distance from
#' infrastructure or disturbance sources in 2 dimensions space.
#'
#' @example examples/zoi_functions_example.R
#'
#' @rdname zoi_functions
#' @export
dist_decay <- function(x, radius = NULL,
                       type = c("exp_decay", "gaussian_decay", "linear_decay",
                                "threshold_decay")[1],
                       zoi_limit = 0.05,
                       origin = 0,
                       oneside = TRUE,
                       ...) {

  if(type %in% c("exp_decay", "exp", "exponential")) {
    return(exp_decay(x = x, radius = radius, zoi_limit = zoi_limit,
                     origin = origin, oneside = oneside, ...))
  }

  if(type %in% c("Gauss", "gauss", "Gaussian", "gaussian", "gaussian_decay",
                 "normal", "Normal", "half_norm", "half_norm_decay")) {
    return(gaussian_decay(x = x, radius = radius, zoi_limit = zoi_limit,
                          origin = origin, oneside = oneside, ...))
  }

  if(type %in% c("Bartlett", "bartlett", "bartlett_decay", "linear",
                 "linear_decay", "tent", "tent_decay")) {
    return(linear_decay(x = x, radius = radius,
                        origin = origin, oneside = oneside, ...))
  }

  if(type %in% c("threshold", "threshold_decay", "step", "step_decay")) {
    return(threshold_decay(x = x, radius = radius,
                           origin = origin, oneside = oneside, ...))
  }
}

#' @rdname zoi_functions
#' @export
threshold_decay <- function(x, radius, intercept = 1, origin = 0, oneside = TRUE) {
  UseMethod("threshold_decay")
}

#' @name zoi_functions
#' @export
threshold_decay.numeric <- function(x, radius, intercept = 1, origin = 0, oneside = TRUE) {
  if(oneside) func <- identity else func <- abs
  ifelse(func(x - origin) < radius, intercept, 0)
}

# possibly have that for RasterLayer as well

#' @name zoi_functions
#' @export
threshold_decay.SpatRaster <- function(x, radius, intercept = 1, origin = 0, oneside = TRUE) {
  if(oneside) func <- identity else func <- abs
  terra::ifel(func(x - origin) < radius, intercept, 0)
}

#' @rdname zoi_functions
#' @export
step_decay <- threshold_decay

#' @rdname zoi_functions
#' @export
bartlett_decay <- function(x, radius, intercept = 1, origin = 0, oneside = TRUE) {
  UseMethod("bartlett_decay")
}

#' @rdname zoi_functions
#' @export
bartlett_decay.numeric <- function(x, radius, intercept = 1, origin = 0, oneside = TRUE) {
  if(oneside) func <- identity else func <- abs
  beta = -intercept/radius
  ifelse(func(x - origin) < radius, intercept + beta * func(x - origin), 0)
}

#' @rdname zoi_functions
#' @export
bartlett_decay.SpatRaster <- function(x, radius, intercept = 1, origin = 0, oneside = TRUE) {
  if(oneside) func <- identity else func <- abs
  beta = -intercept/radius
  terra::ifel(func(x - origin) < radius, intercept + beta * func(x - origin), 0)
}

#' @rdname zoi_functions
#' @export
tent_decay <- bartlett_decay

#' @rdname zoi_functions
#' @export
linear_decay <- bartlett_decay

#' @param sigma `[numeric(1)=NULL]` \cr Standard deviation of the Gaussian decay
#' function. It is related to the Gaussian decay rate \eqn{\lambda} as
#' `lambda = 1/(2*sigma^2)`. Only considered to compute the ZOI
#' for the `gaussian_decay` function when the ZOI radius parameter is null
#' (`radius = NULL`).
#'
#' @param lambda `[numeric(2)=NULL]` \cr For the `gaussian_decay` and `exp_decay`
#' functions, `lambda` is the decay parameter of the Gaussian or exponential decay
#' function. Notice that the interpretation of `lambda` is different depending on the
#' the function -- see details for definitions.
#' For the Gaussian decay function, the value for `lambda` is only considered if both
#' `radius = NULL` and `sigma = NULL`. For the exponential decay function,
#' the value for `lambda` is only considered if both `radius = NULL` and `half_life = NULL`.
#'
#' @rdname zoi_functions
#' @export
gaussian_decay <- function(x, radius = NULL,
                           zoi_limit = 0.05,
                           intercept = 1,
                           lambda = NULL,
                           sigma = NULL,
                           origin = 0, ...) {

  if(!is.null(radius)) {
    lambda = log(1/zoi_limit) / (radius**2)
  } else {
    if(!is.null(sigma)) {
      lambda = 1/(2*sigma**2)
    } else {
      lambda <- lambda
    }

  }

  intercept * exp(- lambda * (x - origin)**2)
}

#' @rdname zoi_functions
#' @export
half_norm_decay <- gaussian_decay

#' @param half_life `[numeric(1)=NULL]` \cr Half life of the exponential decay
#' function, in meters (or map units, for rasters). By definition, the half life is
#' the distance where the exponential decay function reaches 0.5 of its
#' maximum value. For the `exp_decay` function,
#' if the ZOI radius parameter is null (`radius = NULL`), the value of the
#' exponential half life (`half_life = log(2)/lambda`) can be used to parameterize the
#' exponential decay function.
#'
#' @param zoi_hl_ratio `[numeric(1)=NULL]` \cr For the `exp_decay` function,
#' if both the ZOI radius `radius` and `zoi_hl_ratio` are given and
#' `half_life` is `NULL`, this value is used
#' to set the ZOI radius (and `zoi_limit` is ignored).
#' `zoi_hl_ratio` is the ratio between the
#' ZOI radius value and the half life of the exponential function.
#' For instance, if `radius = 1200` and `zoi_hl_ratio = 6`, this means
#' `half_life` is 200. As a consequence, the exponential decay ZOI function
#' decreases to 0.5 at distance 200, and the ZOI radius = 1200
#' is defined as the distance
#' at which the ZOI decreases to 0.5**6 = 0.015625.
#'
#' @rdname zoi_functions
#' @export
exp_decay <- function(x, radius = NULL,
                      zoi_limit = 0.05,
                      intercept = 1,
                      lambda = NULL,
                      origin = 0,
                      oneside = TRUE,
                      half_life = NULL,
                      zoi_hl_ratio = NULL) {

  # define lambda depending on the input parameter
  if(!is.null(radius)) {

    if(is.null(zoi_hl_ratio)) {
      lambda <- log(1/zoi_limit) / radius
    } else {
      half_life <- radius/zoi_hl_ratio
      lambda <- log(2)/half_life
    }

  } else {

    if(!is.null(half_life)) {
      lambda <- log(2)/half_life
    } else {
      lambda <- lambda
    }
  }

  # return function
  if(oneside) func <- identity else func <- abs

  intercept * exp(- lambda * func(x - origin))
}
