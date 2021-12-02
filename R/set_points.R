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
#' If `method = "NLMR"`, the function also uses a raster to define the probabilities of 
#' setting a given point in a certain pixel in space, but this raster is created with
#' a function from the [NLMR] package. The function name is defined by the argument
#' `nlmr_function` and its arguments must be defined as additional parameters to 
#' `set_points()`.
#'
#' TO IMPROVE: implement rasterization with terra package
#'
#' @param n_features `[integer(1)=1000]` \cr Total number of features to spread in space.
#' @param method `[character(1)]{"mobsim", "regular", "random", "raster", "NLMR"}` \cr Method used 
#' to simulate points in space.
#' `mobsim` uses the function [mobsim::sim_thomas_community()] from the [mobsim] package to simulate
#' points. `raster` uses a base raster map as input to define weights and simulate the random points.
#' `NLMR` creates a neutral landscape model using [NLMR] package and uses it as an input base raster. 
#' See `Details` for more information.
#' @param centers `[integer(1)=1]` \cr Number of centers around which the features will be placed.
#' Used only if `method = "mobsim"`.
#' @param width `[numeric(1)=0.05]` \cr Radius of the "patches" of features, around the "patch" centers.
#' Used only if `method = "mobsim"`.
#' @param base_raster `[RasterLayer=NULL]` \cr Base raster to define weights for creating the random points.
#' Used only if `method = "raster"`.
#' @param nlmr_function `[character(1)="nlm_mpd"]` \cr Name of the function from NLMR package used to create
#' the base raster, to be used to define weights for creating the random points.
#' Used only if `method = "NLMR"`.
#' @param res `[numeric(1)=0.1]` \cr Resolution of the output raster.
#' @param extent_x,entent_y `[numeric vector(2)=c(0,1)]` \cr Vector representing the minimum and
#' maximum extent in x and y within which the points should be placed, in the format c(min,max).
#' @param buffer_around `[numeric(1)=0.1]` \cr Size of the buffer around the extent of the landscape,
#' to avoid edge effects when calculating densities using neighborhood analysis.
#' @param return_base_raster `[logical(1)=TRUE]` \cr Whether the base_raster should be returned in
#' the output list. This is `NULL` for `method = "mobsim"`.
#' @param ... Other arguments passed as input to the NLMR functions, defined by the `nlmr_function`
#' argument.
#'
#' @returns A list with three elements: (1) `pts`, the coordinates (x,y) of the simulated points;
#' (2) `rast`, a binary raster containing the landscape, with 1 where there points and NA elsewhere;
#' (3) `base_rast`, the base raster used to weigh the simulation of points. If `method = "mobsim"`,
#' `base_rast` is `NULL`.
#'
#' @example examples/set_points_example.R
#'
#' @export

# function to simulate points in the landscape
set_points <- function(n_features = 1000, 
                       method = c("mobsim", "regular", "random", "raster", "NLMR")[1],
                       centers = 1, width = 0.05,
                       base_raster = NULL,
                       nlmr_function = "nlm_mpd",
                       res = 0.1,
                       extent_x = c(0,1), extent_y = c(0,1),
                       buffer_around = 0,
                       return_base_raster = TRUE,
                       ...) {

  # simulate points according to the method
  if(method == "mobsim") {
    # simulate points with mobsim
    pts <- mobsim::sim_thomas_community(s_pool = 1, n_sim = n_features,
                                        sigma = width, mother_points = centers,
                                        xrange = extent_x, yrange = extent_y)$census[,1:2]
  } else {
    if(method %in% c("NLMR", "raster")) {
      
      # NLMR
      if(method == "NLMR") {
        # simulate points with NLMR
        # get function
        nlm_func <- get(nlmr_function)
        # get nrow and ncol
        ncol = round(abs(diff(extent_x))/res)
        nrow = round(abs(diff(extent_y))/res)
        # simulate landscape
        base_raster <- nlm_func(nrow = nrow, ncol = ncol, resolution = res,
                                ...)
      } 
      
      # simulate points
      pts <- set_points_from_raster(base_raster = base_raster, 
                                    n_features = n_features)
    } else {
      # use set_points_sample
      pts <- set_points_sample(n_features = n_features, type = method,
                               extent_x = extent_x, extent_y = extent_y)
    }
  }

  # extent and res for method "raster"
  if(method == "raster") {
    res = raster::res(base_raster)[1]
    extent_x <- bbox(base_raster)[1,]
    extent_y <- bbox(base_raster)[2,]
  }
  
  # raster
  buff <- buffer_around
  r <- raster::raster(xmn = extent_x[1]-buff, xmx = extent_x[2]+buff,
                      ymn = extent_y[1]-buff, ymx = extent_y[2]+buff,
                      resolution = res)
  # resterize points
  r_pts <- raster::rasterize(pts, r, field = 1)
  
  # return base raster
  if(!return_base_raster) base_raster <- NULL

  # retun
  list(pts = pts, rast = r_pts, base_rast = base_raster)
}
