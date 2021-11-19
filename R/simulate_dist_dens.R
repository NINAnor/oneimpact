#' Simulate landscapes with features and calculate distance and densities of features
#'
#' This function is a combination of [set_points()] and [calc_dist_dens()], to simulate
#' a point pattern representing the distribution of point-type infrastructure in
#' space, rasterize it, and calculate the distance to and density of infrastructure.
#'
#' The function builds upon the function [thomas::sim_thomas_community()] from the [mobsim]
#' package. Originally the function is intended to simulate positions of multiple
#' species in the context of species abundance distribution studies, but it fits
#' well in case of a single species (or point patterns for a single type of feature).
#'
#' TO IMPROVE: implement rasterization with [terra] package
#'
#' @inheritParams set_points
#' @inheritParams calc_dist_dens

#' @returns A RasterBrick with de distance to the nearest feature and the densities for all scales selected.
#'
#' @example examples/simulate_dist_dens_example.R
#'
#' @export

# function to simulate points and calc distance and density rasters
simulate_dist_dens <- function(n_features = 1000, centers = 1,
                               width = 0.05, res = 0.1,
                               extent_x = c(0,1), extent_y = c(0,1),
                               buffer_around = 0,
                               transform_dist = NULL,
                               log_base = exp(1), dist_offset = 1,
                               type_density = c("Gauss", "circle", "rectangle")[1],
                               scale = 100,
                               extent_x_cut = extent_x, extent_y_cut = extent_y) {

  # simulate points
  pts <- set_points(n_features = n_features, centers = centers,
                    width = width, res = res,
                    extent_x = extent_x, extent_y = extent_y,
                    buffer_around = buffer_around)

  # calculate distance and density
  calc_dist_dens(pts$rast, transform_dist = transform_dist,
                 log_base = log_base, dist_offset = dist_offset,
                 type_density = type_density, scale = scale,
                 extent_x_cut = extent_x_cut, extent_y_cut = extent_y_cut)
}
