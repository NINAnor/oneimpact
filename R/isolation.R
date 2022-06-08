#' Isolation and mean isolation in space
#'
#' Measures of isolation and mean isolation to a set of points in
#' space. `isolation` creates random points in a landscape and calculates
#' the nearest neighbor distance from each of them to another set
#' of points passed as input, `x`. `mean_isolation` calculates the
#' average isolation calculated through `isolation`.
#'
#' So far the function only works for a square landscape. In the future
#' we can implement that for polygons or rasters with masks or
#' null cells if necessary, in an approach similar to set_points_sample.
#'
#' @param x `[data.frame]` \cr `data.frame` with (x,y) coordinates in
#' the columns.
#' @param n_rand `[numeric(1)=100]` \cr Number of random points to be
#' created in space, to compute the distance to `x`.
#' @param ext `[numeric(x)=c(0,1)]` Extent of the space within which the
#' random positions should be created c(x or ymin, x or ymax).
#' @param lonlat `[logical(1)=FALSE]` \cr Whether the distance between points
#' should be calculated in an WGS ellipsoid (`lonlat = TRUE`) or on a
#' plane (`lonlat = FALSE`). See [raster::pointDistance] for more details.
#'
#' @return `isolation()` returns the distance from each random point to
#' the nearest neighbor point in `x`. `mean_isolation()` returns the average
#' nearest neighbor distance from all random positions to the points in
#' `x`.
#'
#' @name isolation
#' @export
#'
#' @examples
#' pts <- set_points(n_features = 100, method = "mobsim", centers = 1, width = 0.1)[[1]]
#' isolation(pts)
#' mean_isolation(pts)

# isolation to random points
isolation <- function(x, n_rand = 100, ext = c(0, 1, 0, 1), lonlat = FALSE) {
  # n_rand cannot be the same of the number of points in x
  if(n_rand == nrow(x)) n_rand <- n_rand + 1
  # create random points
  rand <- data.frame(x = runif(n_rand, ext[1], ext[2]), y = runif(n_rand, ext[3], ext[4]))
  # calc dist
  dists  <- raster::pointDistance(x, rand, lonlat = lonlat)
  # min dist (nearest neighbor)
  apply(dists, 2, min)
}

#' @export
#' @rdname isolation
# mean isolation
mean_isolation <- function(x, n_rand = 100, ext = c(0, 1, 0, 1), lonlat = FALSE) {
  mean(isolation(x, n_rand = n_rand, ext = ext, lonlat = lonlat))
}
