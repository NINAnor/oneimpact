#' Calculate density of features at multiple scales
#'
#' This function takes in a raster with locations of infrastructure and calculates
#' a raster (or set of rasters, in case there is more the one value for `scale`)
#' representing the density of features in space (through a spatial filter/neighborhood analysis).
#' The neighborhood analysis is done with the [terra::focal()] function.
#'
#' The neighborhood analysis can be done with different methods. The default is a circular filter
#' (`type = "circle"`), in which case the parameter `scale` corresponds to the radius of the circle
#' centered on the central pixel. Other possibilities are a Gaussian filter
#' (`type = "Gauss"`), in which case scale corresponds to the sigma parameter of a Gaussian
#' functionl and `type = "rectangle"`, in which case the scale corresponds to the
#' size size of the rectangle. For all these methods, these parameters feed the function 
#' [terra::focalMat()] to create the input weight matrix. See [terra::focalMat()] for more
#' details.
#' 
#' If one wants to use their own filter or weight matrix, it is possible to define `type = "mfilter"`
#' and provide a matrix to the parameter `scale` instead, such as one created through the 
#' [terra::focalMat()] or the [oneimpact::create_filter()] functions.
#'
#' TO IMPROVE1: implement with `terra`. WE SHOULD DETECT IF THE INPUT IS RASTER OR TERRA
#'
#' TO IMPROVE2: do the same in communication with GRASS GIS.
#'
#' @param points `[RasterLayer,SpatRaster]` \cr Raster representing locations of features, with 1 where the features
#' are located and NA elsewhere. Can be a [RasterLayer] from [raster] package or a [SpatRaster] from
#' [terra] package.
#' @param type `[character(1)="circle"]{"circle", "Gauss", "rectangle", "mfilter"}` \cr
#' Type of filter used to calculate density. See description for details.
#' @param scale `[numeric(1)=100]` \cr Scale of the neighborhood analysis, used to calculate densities.
#' It can be a single value of a vector of values, in which case several density maps (for each scale)
#' are created. For `type = "circle"`, scale corresponds to the radius of the circle filter. For `type = "Gauss"`,
#' it corresponds to the standard deviation of the Gaussian distribution. If `type = "rectangle"`, it corresponds
#' to the size of the side of a square filter. See [terra::focalMat()] for more details.  
#' If `type = "mfilter"`, scale is not a numeric value but a matrix itself, defined by the user. See description
#' in the details.
#' @param extent_x_cut,entent_y_cut `[numeric vector(2)=c(0,1)]` \cr Vector representing the minimum and
#' maximum extent in x and y for the final output, in the format c(min,max). It is intended to keep only
#' a region of interest but consider the surroundings when calculating densities. The default is to
#' keep the same extent of the input raster.
#' @param na.rm `[logical(1)=FALSE]` \cr Should missing values be removed for filtering calculations?
#' Option for the neighborhood analysis, performed through the [terra::focal()] function.
#' @param plotit `[logical(1)=FALSE]` \cr Should the outputs be plotted along the calculation?
#' @param ... Other arguments to be used within [terra::focal()].
#'
#' @returns A [RasterLayer] or [SpatRaster] (according to the input `points` map) with the density of features. 
#' Alternatively, a `RasterBrick` 
#' or multi-layer `SpatRaster`, if multile `scale` values are given, with the density of features for all scales selected.
#'
#' @example examples/calc_dens_example.R
#'
#' @export
calc_dens <- function(points,
                      type = c("circle", "Gauss", "rectangle", "mfilter")[1],
                      scale = 100,
                      extent_x_cut = terra::ext(points)[c(1,2)],
                      extent_y_cut = terra::ext(points)[c(3,4)],
                      na.rm = FALSE,
                      plotit = FALSE, ...) {

  # check if the input is a terra or raster object
  if(class(points) %in% c("SpatRaster")) {
    use_terra <- TRUE
  } else {
    # we should check if it is raster again here
    use_terra <- FALSE
  }
  
  # check the input is a binary (0,1) map, if it presents only a single value
  r0 <- points
  if(use_terra) {
    if(diff(c(r0@ptr$range_min, r0@ptr$range_max)) == 0) r0 <- terra::classify(r0, cbind(NA, 0)) # binary map
  } else {
    if(diff(c(maxValue(r0), minValue(r0))) == 0) r0 <- raster::reclassify(r0, cbind(NA, 0)) # binary map
  }
  # plot(r0)

  # neighborhood analysis
  if(type == "mfilter") {
    density_r <- terra::focal(r0, w = scale, na.rm = na.rm, ...)
  } else {
    if(length(scale) == 1) {
      gf <- terra::focalMat(r0, d = scale, type = type)
      density_r <- terra::focal(r0, w = gf, na.rm = na.rm, ...)
    } else {
      dens <- sapply(scale, function(x) {
        gf <- terra::focalMat(r0, d = x, type = type)
        terra::focal(r0, w = gf, na.rm = na.rm, ...)
        })
      if(use_terra) density_r <- do.call(c, dens) else
        density_r <- raster::stack(dens)
    }
  }

  # rename
  if(type == "mfilter") names_density = "mfilter" else
    names_density <- paste0("density", scale)
  names(density_r) <- names_density
  if(plotit) plot(density_r)

  # return cropped raster
  if(use_terra)
    terra::crop(density_r, terra::ext(c(extent_x_cut, extent_y_cut)))
  else
    raster::crop(density_r, raster::extent(c(extent_x_cut, extent_y_cut)))
}
