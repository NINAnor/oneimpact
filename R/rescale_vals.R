#' Rescale raster values
#' 
#' @example examples/rescale_vals_example.R
#' 
#' @export
rescale_vals <- function(x, to = c(0, 1), from = NULL) {
  
  # if original miximum and maximum are not given,
  # get them from the layers
  if(is.null(from)) {
    mins <- global(x, "min", na.rm = T)[,1]
    maxs <- global(x, "max", na.rm = T)[,1]
  } else {
    mins <- from[1]
    maxs <- from[2]
  }
  
  # rescale
  x_out <- purrr::pmap(
    list(terra::split(x, 1:nlyr(x)), mins, maxs),
    function(a, b, c) {
      values(a) <- scales::rescale(values(a), to = to, from = c(b, c))
      a
    })
  
  # combine SpatRasters
  x_out <- do.call(c, x_out)
  x_out
}
