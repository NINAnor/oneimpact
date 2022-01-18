#' Calculate the influence from the nearest feature
#'
#' This function takes in a raster with locations of infrastructure and calculates
#' a raster representing the influence to each pixel from the neareast feature of that type
#' of infrastructure. By default, the output measure of influence is the Euclidean distance
#' to the nearest feature. However, this influence measure can be changed to 
#' transformed distance measures (so far, log- and sqrt-distance) or to decay distance functions
#' (so far, exponential decay distance and Bartlett distance). For the decay distance 
#' measures of influence, a zone of influence (zoi) might be used as input for designing 
#' the curve.
#' 
#' The input raster should have positive values in the pixels where infrastructure
#' are located and NA/no-data in all other places. The input raster is supposed to 
#' represent the location of point, line, or polygon infrastructure (e.g. houses, roads, mining areas),
#' but any landscape variable whose representation might be one of those would fit here
#' (e.g. areas of forest or any other habitat type or land cover). We recommend that the
#' input raster is on a metric projection, so that distances and influence measures are
#' based on distance to infrastructure measured in meters.
#' 
#' Explain here better what is the log, sqrt, exp_decay, bartlett, and when a ZOI is defined.
#'
#' TO IMPROVE3: do the same in communication with GRASS GIS.
#' TO IMPROVE4: implement any generic funtion as input to "transform".
#'
#' @param x `[RasterLayer,SpatRaster]` \cr Raster representing locations of features, with value 1 
#' (or any other positive value) where the features are located and NA elsewhere. 
#' Can be a [Raster-class] object (e.g. `RasterLayer`) from [raster-package] or a 
#' [SpatRaster] from [terra-package].
#' The features might be represented by points (e.g. houses, cabins, wind turbines), lines
#' (e.g. roads, power lines, railways), or polygons (e.g. mining areas, urban areas, or
#' any other spatial variable represented as polygons or areas). 
#' 
#' @param zoi `[numeric(1)=NULL]` \cr Zone of Influence (ZoI), in map units meters. 
#' The ZoI is the distance, scale, or buffer size around a feature up to which we consider there is 
#' an effect or influence of an infrastructure or variable. It is considered only when
#' `transform = "bartlett"` or `transform = "exp_decay"`. \cr
#' For the Bartlett influence, it corresponds to the distance beyond which the distance is zero. \cr
#' For the exponential decay influence, the ZoI is used to define the half-life and the lambda 
#' of the exponential decay function, based on the parameter `zoi_hl_ratio`, 
#' that defines the ratio between the ZoI and the half-life. Since the half-life is the value
#' where the exponential decay decreases by `0.5`, a ratio of, for instance, `zoi_hl_ratio = 4` (default)
#' would mean that the ZoI is defined as the value where the exponential decay decreases to `0.5^4 = 0.0625`.
#' In this case, if `zoi = 4000` m, this means that the ZoI is four times higher than the half-life, i.e.
#' `half_life = 1000` and `lambda = log(2)/half_life = 6.93e-4`. The definition of a zone of
#' influence does not imply a cuttoff of the exponential decay function but is only used to define 
#' its parameters, based on the defined `zoi_hl_ratio` parameter.
#' 
#' @param transform `[character(1)=NULL]{"log","sqrt", "exp_decay", "bartlett"}` \cr 
#' By default, NULL, when the measure of influence is the Euclidean distance to the nearest
#' feature. 
#' \itemize{
#'   \item If `log`, the influence measure is the log-distance: `log(euclidean_distance, base = log_base)`.
#'   \item If `sqrt`, the influence measure is the square rooted distance: `sqrt(euclidean_distance)`. 
#'   \item If `exp_decay`, the influence measure is the exponential decay distance: \cr
#' `N_0 * exp(-lambda * euclidean_distance)`. `N_0` and `lambda` are parameters to be defined. 
#' The decayment rate lambda might be defined in terms of the exponential half-life or 
#' a definition of zone of influence (ZoI). 
#'   \item If `bartlett`, the influence measure is a triangular tent-shaped decay distance is returned.
#' }
#' See details below.
#' Other options still to be implemented (such as other functions and a generic user-defined
#' function as input).
#' 
#' @param log_base `[numeric(1)=exp(1)]` \cr Base of the logarithm, if `transform = log`.
#' 
#' @param zoi_hl_ratio `[numeric(1)=6]` \cr Ratio between the zone of influence (ZoI) of the infrastrucure 
#' and the half-life of the exponential decay distance function. It is used to define the 
#' the half-life and decay parameter lambda of the exponential decay function, based on a value of ZoI. 
#' By changing this parameter, one can be more or less flexible in the definition of ZoI for 
#' the exponential decay influence.
#' For example: if `zoi_hl_ratio = 4`, the ZoI is defined as the distance from the nearest infrastructure
#' where the decay distance influence decreases to `0.5^4 = 0.0625`. If `zoi_hl_ratio = 6`, the ZoI
#' is defined as tyhe distance when the exponential decay influence decreases to `0.5^6 = 0.015625`
#' 
#' @param half_life `[numeric(1)=NULL]` \cr Half-life of the exponential decay function, in case
#' `transform = exp_decay`. The exponential decay exponent (lambda) from the exponential function is defined as
#' `lambda = log(2)/half_life`. By definition, when one gets away from the source feature by a 
#' distance interval equals to `half_life`, the magnitude of the exponential decay distance decreases by 1/2. 
#' This means that, for instance, at a distance of `4*half_life` to the nearest feature, the exponential decay 
#' influence has a magnitude of `(1/2)^4 = 1/16 ~ 0.06`. This can be useful to define the 
#' Zone of Influence (ZoI) for exponential decay distances. If the `zoi` parameter is not NULL, the
#' `half_life` parameter is ignored and half-life of the exponential decay is defined by the `zoi` and
#' `zoi_hl_ratio` parameters.
#' 
#' @param exp_decay_parms `[numeric(2)=c(1,0.01)]` \cr Parameters (`N_0`, `lambda`) for the exponential decay 
#' influence, if `transform = exp_decay`. The value of `lambda` defined here is used only if 
#' `zoi = NULL` and `half_life = NULL`, otherwise one of these parameters is used to determine `lambda`.
#' By default, `N_0` is defined as 1, which means the influence is 1 where the infrastructure feature
#' is located, and it decreases as the Euclidean distance from it increases.
#' 
#' @param dist_offset `[numeric(1)=1]` \cr Number to add to the Euclidean distance before transforming it,
#' to avoid `-Inf`/`Inf` values (e.g. in the case of log). It should be a very small value compared to the
#' range of values of Euclidean distance, not to influence any further analyses.
#' 
#' @param extent_x_cut,entent_y_cut `[numeric vector(2)=c(0,1)]` \cr Vectors representing the minimum and
#' maximum extent in x and y for the final output, in the format c(min,max). Used to cut the raster to
#' specific smaller extents. The default is to keep the same extent of the input raster.
#' 
#' @param plotit `[logical(1)=FALSE]` \cr Should the outputs be plotted along the calculation?
#' @param ... \cr Adittional parameters passed to [terra::distance()].
#'
#' @returns A `RasterLayer` (or `SpatRaster`, according to the class of the input object) with the 
#' influence of the nearest feature. By default, it is the Euclidean distance to the nearest feature.
#' Depending on the choice of `transform`, the output distance can be log- or sqrt-transformed distance, 
#' or one can choose to calculate the exponential decay or Bartlett decay influence. 
#' Other types of transformation (e.g. Gaussian?) to be implemented in the future.
#'
#' @example examples/calc_inf_nearest_example.R
#'
#' @export

calc_influence_nearest <- function(
  x,
  zoi = NULL,
  transform = NULL, #c("log", "sqrt", "exp_decay", "bartlett")[1],
  log_base = exp(1),
  zoi_hl_ratio = 4,
  half_life = NULL,
  exp_decay_parms = c(1, 0.01),
  dist_offset = 1,
  extent_x_cut = terra::ext(x)[c(1,2)],
  extent_y_cut = terra::ext(x)[c(3,4)],
  plotit = FALSE, ...) {
  
  # check if the input is a terra or raster object
  if(class(x) %in% c("SpatRaster")) {
    use_terra <- TRUE
  } else {
    if(class(x) %in% c("RasterLayer", "RasterBrick", "RasterStack")) {
      use_terra <- FALSE
    } else {
      classes <- c("SpatRaster", "RasterLayer", "RasterBrick", "RasterStack")
      stop(paste0("Please make sure x is an object of one of these classes: ", 
                  paste(classes, collapse = ","), "."))
    }
  }

  # Euclidean distance
  dist_r <- terra::distance(x, ...)
  
  # transform distance
  if(!is.null(transform))
    if(transform == "log") dist_r <- log(dist_r+dist_offset, base = log_base) else
      if(transform == "sqrt") dist_r <- sqrt(dist_r+dist_offset) else
        if(transform == "exp_decay") {
          # define zoi or half life, depending on which is given as input
          
          # if zoi is given, it is used.
          # if zoi is not given:
          if(is.null(zoi)) { 
            # and half_life is not given either:
            if(is.null(half_life)) {
              # lambda is calculated from exp_decay_parms
              lambda <- exp_decay_parms[2]
            } else {
              # if half_life is given, lambda is calculated from that
              zoi <- half_life * zoi_hl_ratio # not used!
              lambda <- log(2)/half_life
            }
          } else {
            # if zoi is given, it is used to
            # define half_life and lambda
            half_life <- zoi/zoi_hl_ratio
            lambda <- log(2)/half_life
          }
            
          dist_r <- exp_decay_parms[1] * exp(-lambda * dist_r) 
        } else
          if(transform == "bartlett") {
            dist_r <- (1 - (1/zoi)*dist_r)
            zero <- dist_r
            size_landscape <- prod(dim(zero)[c(1:2)])
            values(zero) <- rep(0, size_landscape)
            dist_zero_stack <- c(dist_r, zero)
            
            if(use_terra) {
              dist_r <- terra::app(dist_zero_stack, "max")
            } else {
              dist_zero_stack <- raster::stack(dist_zero_stack)
              dist_r <- raster::calc(dist_zero_stack, max)
            }
            
          } else
            stop("You should select an appropriate transformation method for distance.")

  # rename nearest influence layer
  names(dist_r) <- "influence_nearest"
  # should the result be plotted?
  if(plotit) plot(dist_r)

  # return cropped distance layer
  if(use_terra)
    terra::crop(dist_r, terra::ext(c(extent_x_cut, extent_y_cut)))
  else
    raster::crop(dist_r, raster::extent(c(extent_x_cut, extent_y_cut)))
}
