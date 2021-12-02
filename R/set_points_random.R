#' Simulate random points in 2D
#'
#' This function simulates the coordinates of random points in
#' two dimensions by sampling (x,y) from a uniform distribution
#' between extremes. These extremes (the size of the landscape)
#' are defined by the `extent_x` and `extent_y` parameters.
#' 
#' @param n_features `[integer(1)=1000]` \cr Total number of points to spread in space.
#' @param extent_x,entent_y `[numeric vector(2)=c(0,1)]` \cr Vector representing the minimum and
#' maximum extent in x and y within which the points should be placed, in the format c(min,max).
#' 
#' @return A `data.frame` with the (x,y) random coordinates.
#' 
#' @examples 
#' set_points_random(100)
#' 
set_points_random <- function(n_features = 1000,
                              extent_x = c(0,1), 
                              extent_y = extent_x) {
  # simulate points
  data.frame(x = runif(n_features, extent_x[1], extent_x[2]), y = runif(n_features, extent_y[1], extent_y[2]))
}
