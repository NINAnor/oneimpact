#' Rescale raster values
#'
#' Rescales raster values to a given range (by default, to the range [0,1]).
#'
#' @param x `[SpatRaster]` \cr Raster object or colletion of rasters whose values
#' will be rescaled.
#' @param to `[numeric(2)=c(0,1)]` \cr Range of final values (in the format c(min,max))
#' to which the raster will be rescaled.
#' @param from `[numeric(2)=NULL]` \cr Range of original values (in the format c(min,max))
#' from which the raster will be rescaled. If `NULL`, the minimum and maximum values
#' from `x` are used.
#'
#' @return A raster object (or collection of rasters) with values rescaled to
#' a given range (e.g. to the interval [0,1]).
#'
#' @example examples/raster_rescale_example.R
#'
#' @export
raster_rescale <- function(x, to = c(0, 1), from = NULL) {

  # if original miximum and maximum are not given,
  # get them from the layers
  if(is.null(from)) {
    mins <- terra::global(x, "min", na.rm = T)[,1] # cabins@ptr$range_min?
    maxs <- terra::global(x, "max", na.rm = T)[,1] # cabins@ptr$range_max?
  } else {
    mins <- from[1]
    maxs <- from[2]
  }

  # rescale
  x_out <- mapply(function(a, b, c) {
    terra::values(a) <- scales::rescale(values(a), to = to, from = c(b, c))
    a
  }, a = terra::split(x, 1:terra::nlyr(x)), b = mins, c = maxs)

  # combine SpatRasters
  x_out <- do.call(c, x_out)
  x_out
}
