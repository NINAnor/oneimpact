#' Calculates the influence from the nearest feature and the cumulative influence of infrastructure features
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
#' TO IMPROVE2: do the same in communication with GRASS GIS.
#'
#' @inheritParams calc_influence_nearest
#' @inheritParams calc_influence_cumulative
#'
#' @returns A RasterBrick with de distance to the nearest feature and the densities for all scales selected.
#'
#' @example examples/calc_influence_example.R
#'
#' @export

# function to calculate dist and density
calc_influence <- function(x,
                           zoi,
                           transform_nearest = NULL,
                           type_cumulative = c("circle", "Gauss", "rectangle", "mfilter")[1],
                           extent_x_cut = bbox(x)[1,],
                           extent_y_cut = bbox(x)[2,], 
                           ...) {
  
  # check if the input is a terra or raster object
  if(class(x) %in% c("SpatRaster")) {
    use_terra <- TRUE
  } else {
    if(class(x) %in% c("RasterLayer", "RasterBrick", "RasterStack")) {
      use_terra <- FALSE
    } else {
      classes <- c("SpatRaster", "RasterLayer", "RasterBrick", "RasterStack")
      stop(paste0("Please make sure x is an object of one of these classes: ", 
                  paste(classes, collapse = ","), "."))
    }
  }

  # nearest influence
  nearest_r <- calc_influence_nearest(x = x,
                                      transform = transform_nearest,
                                      extent_x_cut = extent_x_cut,
                                      extent_y_cut = extent_y_cut, 
                                      ...)

  # cumulative influence
  cumulative_r <- calc_influence_cumulative(x = x, zoi = zoi,
                                            type = type_cumulative, 
                                            extent_x_cut = extent_x_cut, 
                                            extent_y_cut = extent_y_cut, 
                                            ...)

  # stack
  if(use_terra) r_stk <- do.call(c, list(nearest_r, cumulative_r)) else
    r_stk <- raster::stack(nearest_r, cumulative_r)
  
  r_stk
}
