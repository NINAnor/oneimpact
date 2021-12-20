#' Calculate distance from the nearest feature
#'
#' This function takes in a raster with locations of infrastructure and calculates
#' a raster representing the distance from each pixel to the neareast feature.
#' The output distance can be transformed (so far, with the options log-distance, sqrt-distance, 
#' exponential decay distance, Bartlett distance).
#'
#' TO IMPROVE1: implement with `terra`.
#'
#' TO IMPROVE2: do the same in communication with GRASS GIS.
#'
#' TO IMPROVE3: Add other possible transformations to distance.
#'
#' @param points `[RasterLayer]` \cr Raster representing locations of features, with 1 where the features
#' are located and NA elsewhere.
#' @param transform_dist `[character(1)=NULL]{"log","sqrt", "exp_decay", "bartlett}` \cr 
#' By default, NULL - distances are no transformed. If `log`, the distances are
#' log-transformed. If `sqrt`, the output is `sqrt(distance)`. If `exp_decay`, the exponential
#' decay distance is calculated. If `bartlett`, a triangular tent-shaped decay distance is returned.
#' See details below.
#' Other options still to be implemented.
#' @param log_base `[numeric(1)=exp(1)]` \cr Base of the logarithm, if `transform_dist = log`.
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
#' @param bartlett `[numeric(1)=NULL]` \cr Zone of Influence (ZoI) of the Bartlett distance, if
#' `transform_dist = bartlett`. It corresponds to the distance beyonf which the distance is zero.
#' @param dist_offset `[numeric(1)=1]` \cr Number to add to distance before transforming it,
#' to avoid `-Inf`/`Inf` values (e.g. in the case of log). It should be a very small value compared to the
#' range of values of distance.
#' @param extent_x_cut,entent_y_cut `[numeric vector(2)=c(0,1)]` \cr Vectors representing the minimum and
#' maximum extent in x and y for the final output, in the format c(min,max). The default is to
#' keep the same extent of the input raster.
#' @param plotit `[logical(1)=FALSE]` \cr Should the outputs be plotted along the calculation?
#'
#' @returns A RasterLayer with the distance to the nearest feature. Depending on the choice of
#' `transform_dist`, the output distance can be log- or sqrt-transformed, or one can choose to calculate
#' the exponential decay or Bartlett decay distance. Other types of transformation
#' to be implemented in the future.
#'
#' @example examples/calc_dist_example.R
#'
#' @export

calc_dist <- function(points,
                      transform_dist = NULL, #c("log", "sqrt", "exp_decay", "bartlett")[1],
                      log_base = exp(1),
                      half_life = NULL, #log(2)/0.01,
                      exp_decay_parms = c(1, 0.01),
                      bartlett_zoi = NULL,
                      dist_offset = 1,
                      extent_x_cut = bbox(points)[1,],
                      extent_y_cut = bbox(points)[2,],
                      use_terra = TRUE,
                      plotit = FALSE) {

  # distance
  dist_r <- terra::distance(points)
  
  # transform distance
  if(!is.null(transform_dist))
    if(transform_dist == "log") dist_r <- log(dist_r+dist_offset, base = log_base) else
      if(transform_dist == "sqrt") dist_r <- sqrt(dist_r+dist_offset) else
        if(transform_dist == "exp_decay") {
          if(is.null(half_life)) lambda <- exp_decay_parms[2] else lambda <- log(2)/half_life
          dist_r <- exp_decay_parms[1] * exp(-lambda * dist_r) 
        } else
          if(transform_dist == "bartlett") {
            dist_r <- (1 - (1/bartlett_zoi)*dist_r)
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

  names(dist_r) <- "distance"
  if(plotit) plot(dist_r)

  terra::crop(dist_r, extent(c(extent_y_cut, extent_x_cut)))
}
