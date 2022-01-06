#' Calculate distance from the nearest feature and density of features
#'
#' This function takes in a raster with locations of infrastructure and calculates (1)
#' a raster representing the distance from each pixel to the neareast feature and (2)
#' a raster (or set of rasters, in case there is more the one value for `scale`)
#' representing the density of features in space (through a spatial filter/neighborhood analysis).
#' The neighborhood analysis is done with the [raster::focal()] function.
#'
#' The neighborhood analysis can be done with different methods. The default is a Gaussian filter
#' (`type_density = "Gauss"`), in which case scale corresponds to the sigma paramater of the Gaussian
#' filter. If `type_density = "circle"` or `type_density = "rectangle"`, the scale corresponds to the
#' radius of the circle or width of the rectangle, respectively. See [raster::focalWeight()] for more
#' details.
#'
#' TO IMPROVE1: implement with `terra`.
#'
#' TO IMPROVE2: do the same in communication with GRASS GIS.
#'
#' TO IMPROVE3: Add other possible transformations to distance.
#'
#' @inheritParams calc_dist
#' @inheritParams calc_dens
#'
#' @returns A RasterBrick with de distance to the nearest feature and the densities for all scales selected.
#'
#' @example examples/calc_dist_dens_example.R
#'
#' @export

# function to calculate dist and density
calc_dist_dens <- function(points,
                           transform_dist = NULL,
                           log_base = exp(1), dist_offset = 1,
                           half_life = NULL,
                           exp_decay_parms = c(1, 0.01),
                           bartlett_zoi = NULL,
                           type_density = c("Gauss", "circle", "rectangle")[1],
                           scale = 100,
                           extent_x_cut = bbox(points)[1,],
                           extent_y_cut = bbox(points)[2,],
                           plotit = FALSE, ...) {

  # distance
  dist_r <- calc_dist(points = points, transform_dist = transform_dist, 
                      log_base = log_base, half_life = half_life, 
                      exp_decay_parms = exp_decay_parms, 
                      bartlett_zoi = bartlett_zoi,
                      dist_offset = dist_offset,
                      extent_x_cut = extent_x_cut, 
                      extent_y_cut = extent_y_cut, 
                      plotit = plotit)

  # density
  density_r <- calc_dens(points, type_density, scale,
                         extent_x_cut, extent_y_cut,
                         plotit, ...)

  r_stk <- raster::stack(dist_r, density_r)
  r_stk
}
