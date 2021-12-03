#' Calculate density of features at multiple scales
#'
#' This function takes in a raster with locations of infrastructure and calculates
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
#' TO IMPROVE1: implement with `terra`.
#'
#' TO IMPROVE2: do the same in communication with GRASS GIS.
#'
#' @param points `[RasterLayer]` \cr Raster representing locations of features, with 1 where the features
#' are localted and NA elsewhere.
#' @param type_density `[character(1)="Gauss"]{"Gauss", "circle", "rectangle"}` \cr
#' Type of filter used to calculate density. See description for details.
#' @param scale `[numeric(1)=100]` \cr Scale of the neighborhood analysis, used to calculate densities.
#' It can be a single value of a vector of values, in which case several density maps (for each scale)
#' are created.
#' @param extent_x_cut,entent_y_cut `[numeric vector(2)=c(0,1)]` \cr Vector representing the minimum and
#' maximum extent in x and y for the final output, in the format c(min,max). It is intended to keep only
#' a region of interest but consider the surroundings when calculating densities. The default is to
#' keep the same extent of the input raster.
#' @param `plotit [logical(1)=FALSE]` \cr Should the outputs be plotted along the calculation?
#'
#' @returns A RasterBrick with the density of features for all scales selected.
#'
#' @example examples/calc_dens_example.R
#'
#' @export
calc_dens <- function(points,
                      type_density = c("Gauss", "circle", "rectangle")[1],
                      scale = 100,
                      extent_x_cut = bbox(points)[1,],
                      extent_y_cut = bbox(points)[2,],
                      plotit = FALSE, ...) {

  # density
  r0 <- points
  if(diff(c(maxValue(r0), minValue(r0))) == 0) r0[is.na(r0[])] <- 0 # binary map
  # plot(r0)

  # Gaussian weight
  if(length(scale) == 1) {
    gf <- raster::focalWeight(r0, d = scale, type = type_density)
    density_r <- raster::focal(r0, w = gf, ...)
  } else {
    dens <- sapply(scale, function(x) {
      gf <- raster::focalWeight(r0, d = x, type = type_density)
      raster::focal(r0, w = gf, ...)
    })
    density_r <- raster::stack(dens)
  }

  names_density <- paste0("density", scale)
  names(density_r) <- names_density
  if(plotit) plot(density_r)

  raster::crop(density_r, extent(c(extent_y_cut, extent_x_cut)))
}
