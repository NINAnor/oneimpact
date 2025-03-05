#' Simulate points in a landscape
#'
#' This function simulates point patterns in space and rasterize them.
#' The idea is to mimic the spatial distribution of point-type infrastructure,
#' such as houses, cabins, or turbines, for instance.
#' The function returns a list with the position of the points and a binary raster
#' with 1 where there are points and NA elsewhere. If created with a raster to define
#' the weights, this base raster is also returned in the output.
#'
#' If `method = "mobsim"`, the function builds upon the function
#' [mobsim::sim_thomas_community()] from the [mobsim]
#' package. Originally the function is intended to simulate positions of multiple
#' species in the context of species abundance distribution studies, but it fits
#' well in case of a single species (or point patterns for a single type of feature).
#' In this case, the points are simulated based on the number of centers/patches of points
#' and their width.
#'
#' If `method = "raster"`, the function uses an input raster (defined by the argument
#' `base_raster`) to define the probabilities of setting a given point in a certain pixel
#' in space.
#'
#' @param n_features `[integer(1)=1000]` \cr Total number of features to spread in space.
#' @param method `[character(1)]{"mobsim", "regular", "random", "raster"}` \cr Method used
#' to simulate points in space.
#' `mobsim` uses the function [mobsim::sim_thomas_community()] from the [mobsim] package to simulate
#' points. `raster` uses a base raster map as input to define weights and simulate the random points.
#' @param centers `[integer(1)=1]` \cr Number of centers around which the features will be placed.
#' Used only if `method = "mobsim"`.
#' @param width `[numeric(1)=0.05]` \cr Mean distance between each of the features in a cluster
#' and the center of the cluster. Used only if `method = "mobsim"`.
#' @param base_raster `[RasterLayer=NULL]` \cr Base raster to define weights for creating the random points.
#' Used only if `method = "raster"`.
#' @param point_coordinates `[data.frame=NULL]` \cr `data.frame` with (x,y) columns with coordinates
#' already taken from elsewhere. This option is intended for when the points' coordinates were already
#' generated or taken from a real landscape. In this case, no points are simulated and they are
#' just rasterized (so that distances or other derived variables might be calculated).
#' @param res `[numeric(1)=0.1]` \cr Resolution of the output raster.
#' @param extent_x,entent_y `[numeric vector(2)=c(0,1)]` \cr Vector representing the minimum and
#' maximum extent in x and y within which the points should be placed, in the format c(min,max).
#' @param buffer_around `[numeric(1)=0.1]` \cr Size of the buffer around the extent of the landscape,
#' to avoid edge effects when calculating densities using neighborhood analysis.
#' @param return_base_raster `[logical(1)=TRUE]` \cr Whether the base_raster should be returned in
#' the output list. This is `NULL` for `method = "mobsim"`.
#' @param use_terra `[logical(1)=TRUE]` \cr If `TRUE` (default), the `rast` element
#' created from the points is a `SpatRaster` object
#' from `terra` package is created. If `FALSE`, it is a `RasterLayer` from `raster` package
#' is created.
#' @param crs `[character(1)]` \cr Specification for the coordinate reference system
#' of the `rast` object created from the points. Default is
#' `"+proj=utm +zone=1 +datum=WGS84"`.
#' argument.
#'
#' @returns A list with three elements: (1) `pts`, the coordinates (x,y) of the simulated points;
#' (2) `rast`, a binary raster containing the landscape, with 1 where there points and NA elsewhere;
#' (3) `base_rast`, the base raster used to weigh the simulation of points. If `method = "mobsim"`
#' or `"regular"` or `"random"`, `base_rast` is `NULL`.
#'
#' @example examples/set_points_example.R
#'
#' @export

# function to simulate points in the landscape
set_points <- function(n_features = 1000,
                       method = c("mobsim", "regular", "random", "raster")[1],
                       centers = 1, width = 0.05,
                       base_raster = NULL,
                       point_coordinates = NULL,
                       res = 0.1,
                       extent_x = c(0,1), extent_y = c(0,1),
                       buffer_around = 0,
                       return_base_raster = TRUE,
                       use_terra = TRUE,
                       crs = "") {

  # get point coordinates if they were taken from elsewhere
  if(!is.null(point_coordinates)) {
    pts <- point_coordinates
  } else {
    # simulate points according to the method
    if(method == "mobsim") {
      # simulate points with mobsim
      pts <- mobsim::sim_thomas_community(s_pool = 1, n_sim = n_features,
                                          sigma = width, mother_points = centers,
                                          xrange = extent_x, yrange = extent_y)$census[,1:2]
    } else {
      if(method %in% c("raster")) {

        # simulate points
        pts <- oneimpact::set_points_from_raster(base_raster = base_raster,
                                                 n_features = n_features)
      } else {
        # use set_points_sample
        pts <- oneimpact::set_points_sample(n_features = n_features, type = method,
                                            extent_x = extent_x, extent_y = extent_y)
      }
    }
  }

  # extent and res for method "raster"
  if(method == "raster") {
    res <- terra::res(base_raster)[1]
    extent_x <- terra::ext(base_raster)[1:2]
    extent_y <- terra::ext(base_raster)[3:4]
  }

  # raster
  buff <- buffer_around
  if(use_terra) {
    if(crs == "") crs <- "+proj=utm +zone=1 +datum=WGS84" # example of planar crs
    r <- terra::rast(xmin = extent_x[1]-buff, xmax = extent_x[2]+buff,
                     ymin = extent_y[1]-buff, ymax = extent_y[2]+buff,
                     resolution = res, crs = crs)
    # resterize points
    r_pts <- as.matrix(pts) |>
      terra::vect() |>
      terra::rasterize(r, field = 1)
  } else {
    r <- raster::raster(xmn = extent_x[1]-buff, xmx = extent_x[2]+buff,
                        ymn = extent_y[1]-buff, ymx = extent_y[2]+buff,
                        resolution = res)
    # resterize points
    r_pts <- raster::rasterize(pts, r, field = 1)
  }

  # return base raster
  if(!return_base_raster) base_raster <- NULL

  # retun
  list(pts = pts, rast = r_pts, base_rast = base_raster)
}
