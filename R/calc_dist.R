#' Calculate distance from the nearest feature
#'
#' This function takes in a raster with locations of infrastructure and calculates
#' a raster representing the distance from each pixel to the neareast feature.
#' The output distance can be transformed (log, sqrt).
#'
#' TO IMPROVE1: implement with `terra`.
#'
#' TO IMPROVE2: do the same in communication with GRASS GIS.
#'
#' TO IMPROVE3: Add other possible transformations to distance.
#'
#' @param points `[RasterLayer]` \cr Raster representing locations of features, with 1 where the features
#' are located and NA elsewhere.
#' @param transform_dist `[character(1)=NULL]{"log","sqrt"}` \cr By default, NULL. If "log", the distances are
#' log-transformed. If "sqrt", the output is `sqrt(distance)`. Other options still to be implemented.
#' @param log_base `[numeric(1)=exp(1)]` \cr Base of the logarithm, if `log_base` is TRUE.
#' @param dist_offset `[numeric(1)=1]` \cr Number to add to distance before transforming it,
#' to avoid `-Inf`/`Inf` values (e.g. in the case of log). It should be a very small value compared to the
#' range of values of distance.
#' @param extent_x_cut,entent_y_cut `[numeric vector(2)=c(0,1)]` \cr Vectors representing the minimum and
#' maximum extent in x and y for the final output, in the format c(min,max). The default is to
#' keep the same extent of the input raster.
#' @param plotit `[logical(1)=FALSE]` \cr Should the outputs be plotted along the calculation?
#'
#' @returns A RasterLayer with the distance to the nearest feature. Depending on the choice of
#' `transform_dist`, the output distance can be log- or sqrt-transformed. Other types of transformation
#' to be implemented.
#'
#' @example examples/calc_dist_example.R
#'
#' @export

calc_dist <- function(points,
                      transform_dist = NULL, #c(NA, "log", "sqrt", "exp_decay")[1],
                      log_base = exp(1),
                      exp_hl = NULL, #log(2)/0.01,
                      exp_decay_parms = c(1, 0.01),
                      dist_offset = 1,
                      extent_x_cut = bbox(points)[1,],
                      extent_y_cut = bbox(points)[2,],
                      plotit = FALSE) {

  # distance
  dist_r <- raster::distance(points)
  if(!is.null(transform_dist))
    if(transform_dist == "log") dist_r <- log(dist_r+dist_offset, base = log_base) else
      if(transform_dist == "sqrt") dist_r <- sqrt(dist_r+dist_offset) else
        if(transform_dist == "exp_decay") {
          if(is.null(exp_hl)) alfa <- exp_decay_parms[2] else alfa <- log(2)/exp_hl
          dist_r <- exp_decay_parms[1] * exp(-alfa * dist_r) 
        } else
            stop("You should select an appropriate transformation method for distance.")

  names(dist_r) <- "distance"
  if(plotit) plot(dist_r)

  raster::crop(dist_r, extent(c(extent_y_cut, extent_x_cut)))
}
