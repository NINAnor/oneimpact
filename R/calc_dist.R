#' Calculate distance from the nearest feature
#'
#' This function takes in a raster with locations of infrastructure and calculates
#' a raster representing the distance from each pixel to the neareast feature.
#' The output distance can be transformed (so far, with the options log-distance, sqrt-distance, 
#' exponential decay distance, Bartlett distance).
#' 
#' Explain here better what is the log, sqrt, exp_decay, bartlett, and when a ZOI is defined.
#'
#' TO IMPROVE1: implement with `terra`. Seems Ok. WE SHOULD DETECT IF THE INPUT IS RASTER OR TERRA
#' 
#' TO IMPROVE2: define ZOI parameters, defined for exp decay given a certain ratio.
#'
#' TO IMPROVE3: do the same in communication with GRASS GIS.
#'
#' @param x `[RasterLayer,SpatRaster]` \cr Raster representing locations of features, with 1 where the features
#' are located and NA elsewhere. Can be a [RasterLayer] from [raster] package or a [SpatRaster] from
#' [terra] package.
#' @param transform_dist `[character(1)=NULL]{"log","sqrt", "exp_decay", "bartlett}` \cr 
#' By default, NULL - distances are no transformed. If `log`, the distances are
#' log-transformed. If `sqrt`, the output is `sqrt(distance)`. If `exp_decay`, the exponential
#' decay distance is calculated. If `bartlett`, a triangular tent-shaped decay distance is returned.
#' See details below.
#' Other options still to be implemented.
#' @param zoi `[numeric(1)=NULL]` \cr Zone of Influence (ZoI), in meters. The ZoI is distance or scale up to 
#' which we consider there is effect of an infrastructure or variable. It is considered only when
#' `transform_dist = "bartlett"` or `transform_dist = "exp_decay"`. 
#' For the Bartlett distance, it corresponds to the distance beyond which the distance is zero.
#' For the exponential decay distance, the ZoI is used to define the half_life of the exponential decay
#' function, based on the parameter `zoi_hl_ratio`, that defines the ratio between the ZoI and the half life.
#' For instance, if the `zoi = 4000` m and `zoi_hl_ratio = 4`, the ZoI corresponds to the distance at which
#' the exponential decay decreases to `0.5**4 = 0.0625`. In such a case, the half life is defined as
#' `half_life = 1000`.
#' @param log_base `[numeric(1)=exp(1)]` \cr Base of the logarithm, if `transform_dist = log`.
#' @param zoi_hl_ratio `[numeric(1)=6]` \cr Ratio between the ZoI and the half life of the exponential decay
#' distance function. It is used to define the ZoI for the exponential decay function. For instance, if 
#' `half_life = 1000` and `zoi_hl_ratio = 4`, the ZoI will be 4000 m (when the exponential decay decrease to 
#' `0.5**4 = 0.0625`.
#' @param half_life `[numeric(1)=NULL]` \cr Half life of the exponential decay function, in case
#' `transform_dist = exp_decay`. The lambda exponent from the exponential funcion is defined as
#' `lambda = log(2)/half_life`. By definition, t each distance interval equals to `half_life` from the 
#' source features, the magnitude of the exponential decay distance decreases by 1/2. This means that,
#' for instance, at a distance of `4*half_life` to the nearest feature, the exponential decay distance
#' has a magnitude of 1/16 ~ 0.06. This can be useful to define the Zone of Influence for exponential
#' decay distances.
#' @param exp_decay_parms `[numeric(2)=c(1,0.01)]` \cr Parameters (`N_0`, `lambda`) for the exponential decay 
#' distance, if `transform_dist = exp_decay`. The value of `lambda` define here is used only if `half_life = NULL`,
#' otherwise the value of `half_life` is used to determine `lambda`.
#' @param dist_offset `[numeric(1)=1]` \cr Number to add to distance before transforming it,
#' to avoid `-Inf`/`Inf` values (e.g. in the case of log). It should be a very small value compared to the
#' range of values of distance.
#' @param extent_x_cut,entent_y_cut `[numeric vector(2)=c(0,1)]` \cr Vectors representing the minimum and
#' maximum extent in x and y for the final output, in the format c(min,max). The default is to
#' keep the same extent of the input raster.
#' @param plotit `[logical(1)=FALSE]` \cr Should the outputs be plotted along the calculation?
#'
#' @returns A `RasterLayer` (or `SpatRaster`, depending on the input) with the distance to the nearest feature. 
#' Depending on the choice of `transform_dist`, the output distance can be log- or sqrt-transformed, or one can choose to calculate
#' the exponential decay or Bartlett decay distance. Other types of transformation (e.g. Gaussian?)
#' to be implemented in the future.
#'
#' @example examples/calc_dist_example.R
#'
#' @export

calc_dist <- function(x,
                      transform_dist = NULL, #c("log", "sqrt", "exp_decay", "bartlett")[1],
                      zoi = NULL,
                      log_base = exp(1),
                      zoi_hl_ratio = 4,
                      half_life = NULL, #log(2)/lambda,
                      exp_decay_parms = c(1, 0.01),
                      dist_offset = 1,
                      extent_x_cut = terra::ext(x)[c(1,2)],
                      extent_y_cut = terra::ext(x)[c(3,4)],
                      plotit = FALSE) {
  
  # check if the input is a terra or raster object
  if(class(x) %in% c("SpatRaster")) {
    use_terra <- TRUE
  } else {
    # we should check if it is raster again here
    use_terra <- FALSE
  }

  # distance
  dist_r <- terra::distance(x)
  
  # transform distance
  if(!is.null(transform_dist))
    if(transform_dist == "log") dist_r <- log(dist_r+dist_offset, base = log_base) else
      if(transform_dist == "sqrt") dist_r <- sqrt(dist_r+dist_offset) else
        if(transform_dist == "exp_decay") {
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
          if(transform_dist == "bartlett") {
            dist_r <- (1 - (1/zoi)*dist_r)
            zero <- dist_r
            values(zero) <- 0
            dist_zero_stack <- c(dist_r, zero)
            
            if(use_terra) {
              dist_r <- terra::app(dist_zero_stack, "max")
            } else {
              dist_zero_stack <- raster::stack(dist_zero_stack)
              dist_r <- raster::calc(dist_zero_stack, max)
            }
            
          } else
            stop("You should select an appropriate transformation method for distance.")

  # rename distance layer
  names(dist_r) <- "distance"
  if(plotit) plot(dist_r)

  # return cropped distance layer
  if(use_terra)
    terra::crop(dist_r, terra::ext(c(extent_x_cut, extent_y_cut)))
  else
    raster::crop(dist_r, raster::extent(c(extent_x_cut, extent_y_cut)))
}
