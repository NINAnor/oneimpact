#' Simulate points in a landscape
#'
#' This function simulate points patterns in space and rasterize them.
#' The idea is to mimic the spatial distribution of point-type infrastructure,
#' such as houses, cabins, or turbines, for instance.
#' The function returns a list with the position of the points and a binary raster
#' with 1 where there are points and NA elsewhere.
#'
#' The function builds upon the function [mobsim::sim_thomas_community()] from the [mobsim]
#' package. Originally the function is intended to simulate positions of multiple
#' species in the context of species abundance distribution studies, but it fits
#' well in case of a single species (or point patterns for a single type of feature).
#'
#' TO IMPROVE: implement rasterization with terra package
#'
#' @param n_features `[integer(1)=1000]` \cr Total number of features to spread in space.
#' @param centers `[integer(1)=1]` \cr Number of centers around which the features will be placed.
#' @param width `[numeric(1)=0.05]` \cr Radius of the "patches" of features, around the "patch" centers
#' @param res `[numeric(1)=0.1]` \cr Resolution of the output raster.
#' @param extent_x,entent_y `[numeric vector(2)=c(0,1)]` \cr Vector representing the minimum and
#' maximum extent in x and y within which the points should be placed, in the format c(min,max).
#' @param buffer_around `[numeric(1)=0.1]` \cr Size of the buffer around the extent of the landscape,
#' to avoid edge effects when calculating densities using neighborhood analysis.
#'
#' @returns A list with two elements: (1) `pts`, the coordinates (x,y) of the simulated points;
#' (2) `rast`, a binary raster containing the landscape, with 1 where there points and NA elsewhere.
#'
#' @example examples/set_points_example.R
#'
#' @export

# function to simulate points in the landscape
set_points <- function(n_features = 1000, centers = 1,
                       width = 0.05, res = 0.1,
                       extent_x = c(0,1), extent_y = c(0,1),
                       buffer_around = 0) {

  # simulate points
  pts <- mobsim::sim_thomas_community(s_pool = 1, n_sim = n_features,
                                      sigma = width, mother_points = centers,
                                      xrange = extent_x, yrange = extent_y)$census[,1:2]

  # raster
  buff <- buffer_around
  r <- raster::raster(xmn = extent_x[1]-buff, xmx = extent_x[2]+buff,
                      ymn = extent_y[1]-buff, ymx = extent_y[2]+buff,
                      resolution = res)
  # resterize points
  r_pts <- raster::rasterize(pts, r, field = 1)

  # retun
  list(pts = pts, rast = r_pts)
}
