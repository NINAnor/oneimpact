#' Calculates the zone of influence from the nearest feature
#' and the cumulative zone of influence of multiple features
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
#' @inheritParams calc_zoi_nearest
#' @inheritParams calc_zoi_cumulative
#' @param zoi_metric `[character(1)="all"]{"all", "nearest", "cumulative"}` \cr
#' Which metric of zone of influence should be computed. Either `"all"`, `"nearest"`,
#' or `"cumulative"`.
#'
#' @returns A RasterBrick with de distance to the nearest feature and the densities for all scales selected.
#'
#' @seealso [oneimpact::calc_zoi_nearest()], [oneimpact::calc_zoi_cumulative()]
#'
#' @example examples/calc_zoi_example.R

# function to calculate dist and density
calc_zoi <- function(x,
                     zoi_metric = c("all", "nearest", "cumulative")[1],
                     where = c("R", "GRASS")[1],
                     ...) {

  # check if the input is a terra or raster object
  if(where == "R")
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
  if(zoi_metric %in% c("all", "near", "nearest", "dist", "distance")) {
    nearest_r <- calc_zoi_nearest(x = x,
                                  ...)
  }

  # cumulative influence
  if(zoi_metric %in% c("all", "cum", "cumulative", "density")) {
    cumulative_r <- calc_zoi_cumulative(x = x,
                                        ...)
  }

  # stack
  if(zoi_metric == "all") {
    if(where == "R") {
      if(use_terra) r_stk <- do.call(c, list(nearest_r, cumulative_r)) else
        r_stk <- raster::stack(nearest_r, cumulative_r)
    } else {
      r_stk <- do.call(c, list(nearest_r, cumulative_r))
    }
  } else {
    if(zoi_metric %in% c("near", "nearest", "dist", "distance")) {
      r_stk <- nearest_r
    } else {
      r_stk <- cumulative_r
    }
  }

  r_stk
}
