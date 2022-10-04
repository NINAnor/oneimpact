#' Create filters or kernel matrices for raster neighborhood analyses
#'
#' This function creates matrices of weights following different
#' functions to be used in neighborhood analyses for rasters. In the context of
#' cumulative impact analysis, they represent the Zone of Influence (ZoI) of each
#' infrastructure point/pixel, to be used to calculate the cumulative ZoI.
#' It is possible to export these matrices as text files, for use with external
#' software such as the `r.mfilter` module within GRASS GIS.
#'
#' The function creates \eqn{n} x \eqn{n} ZoI or weight matrices based on
#' functions with different shapes and parameterized with the ZoI radius, where
#' \eqn{n} is the dimension of the matrix.
#' For some functions (e.g. threshold decay, linear decay),
#' the size of the matrix is defined by the ZoI radius, in meters,
#' given the intended resolution (parameter `r`), potentially adding new lines
#' and columns with value zero to keep \eqn{n} an odd number.
#' For non-vanishing function (e.g. exponential or Gaussian decay),
#' even though the function is parameterized with the ZoI radius the size of
#' the matrix can go beyond this radius. In this case, the size of the matrix
#' \eqn{n} is defined either by a minimum intensity function value
#' (parameter `min_intensity`) or by a maximum distance for
#' the matrix radius (parameter `min_dist`, which can be set to be the `radius`).
#' Keeping \eqn{n} at a reasonable size guarantees that the neighborhood
#' analysis using such input weight matrices is computationally feasible.
#'
#' Possible future implementation: compare results with
#' [smoothie::kernel2dsmooth()] and [smoothie::kernel2dmeitsjer()],
#' maybe wrap some options here.
#'
#' @param r `[numeric,SpatRaster,RasterLayer]` \cr Either a numeric value
#' corresponding to the resolution (pixel size) that each pixel in the filter matrix
#' should correspond to; or a raster object (`SpatRaster` from the [terra]
#' package or `RasterLayer`, `RasterBrick`, or `RasterStack` from the
#' `raster` package) from which such resolution can be extracted.
#'
#' @param radius `[numeric(1)=NULL]` \cr Zone of Influence (ZoI) radius,
#' in map units (preferentially meters).
#' The ZoI radius is the distance, scale, or buffer size around a
#' feature up to which we consider there is
#' an effect or influence of an infrastructure or variable. In `filter_create`,
#' the interpretation of the
#' radius differ depending on the shape of the zoi (parameter `type`):
#' - For the circle neighborhood (`type = "circle"` or `type = "threshold"`
#' or `type = "step"`), the `radius` corresponds to the radius
#' (or threshold) of the circle, beyond which the filter is zero.
#' - For the rectangular neighborhood (`type = "rectangle"` or `type = "box"`),
#' the `radius` corresponds to half the size of the square size, or
#' `square size = 2*radius`. For a rectangular filter with different size
#' of the sides, use [terra::focal()] (but
#' please note the interpretation of the parameters is different).
#' - For the Bartlett neighborhood (`type = "bartlett"` or
#' `type = "linear_decay"` or `type = "tent_decay"`),
#' the `radius` corresponds to the distance beyond which the filter is zero.
#' - For the exponential decay neighborhood (`type = "exp_decay"`) and the
#' Gaussian decay neighborhood (`type = "Gauss"` or `type = "gaussian_decay"`),
#' the `radius` corresponds to the distance where the exponential decay
#' function goes below a given limit distance defined by
#' `zoi_limit`. See [oneimpact::zoi_functions()] for details.
#' - If `radius = NULL`, the exponential or gaussian decay matrices are
#' defined based on other parameters -- see below. This option will raise an
#' error for the other types of filters.
#'
#' @param zoi_limit `[numeric(1)=0.05]` \cr For non-vanishing filters
#' (e.g. `exp_decay`, `gaussian_decay`), this value is used to set the relationship
#' between the ZoI radius and the decay functions:
#' `radius` is defined as the minimum distance `x` at which the ZoI assumes values
#' below `zoi_limit`. The default is 0.05. This parameter is used only
#' if `radius` is not `NULL`.
#'
#' @param type `[character(1)="exp_decay"]{"exp_decay", "bartlett", "circle",
#' "threshold_decay", "gaussian_decay", "Gauss", "rectangle"}` \cr
#' Shape of the Zone of Influence of weight matrix. It can be any of:
#' - `"circle"`, `"threshold"`, `"threshold_decay"`, `"step"` or `"step_decay"`
#' for a threshold decay ZoI;
#' - `"exp_decay"` for exponential decay ZoI;
#' - `"Gauss"`, `"gaussian"`, or `"gaussian_decay"` for Gaussian decay ZoI;
#' - `"bartlett"`, `"bartlett_decay"`, `"linear_decay"`, or `"tent_decay"`
#' for linear decay ZoI;
#' - `"rectangle"` or `"box"` for a rectangular ZoI.
#' There might be some correspondence between the weight matrix `type`
#' in `filter_create` and other similar functions (e.g. `type = "rectangle"`
#' and `type = "boxcar"` in [smoothie::kernel2dmeitsjer()] or
#' `type = "Gauss"` in [terra::focalMat()] with parameter
#' `type = "gauss"` n [smoothie::kernel2dmeitsjer]()); however, the
#' interpretation of the parameters used to
#' define these matrices is different between functions.
#'
#' @param half_life `[numeric(1)=NULL]` \cr Half life of the exponential decay
#'  function, in meters. By definition, the half life is
#'  the distance where the exponential decay function reaches 0.5 of its
#'  maximum value. For the `exp_decay` function,
#'  if the ZoI radius parameter is null (`radius = NULL`), the value of the
#'  exponential half life (`half_life = log(2)/lambda`) can used to parameterize the
#'  exponential decay function. See details in [oneimpact::zoi_functions()].
#' @param zoi_hl_ratio `[numeric(1)=6]` \cr For the `exp_decay` function,
#' if both the ZoI radius `radius` and `zoi_hl_ratio` are given and
#' `half_life` is `NULL`, this value is used
#' to set the ZoI radius (and `zoi_limit` is ignored).
#' `zoi_hl_ratio` is the ratio between the
#' ZoI radius value and the half life of the exponential function.
#' For instance, if `radius = 1200` and `zoi_hl_ratio = 6`, this means
#' `half_life` is 200. As a consequence, the exponential decay ZoI function
#' decreases to 0.5 at distance 200, and the ZoI radius = 1200
#' is defined as the distance
#' at which the ZoI decreases to 0.5**6 = 0.015625.
#' @param min_intensity `[numeric(1)=0.01]` \cr Minimum intensity of the
#' exponential and Gaussian decay functions to
#' define the radius of the window that define the filter.
#' @param max_dist `[numeric(1)=50000]` \cr Maximum size (in meters) to
#' define the radius of the window that defines the filter. Only
#' applicable for exponential and Gaussian decay functions.
#' @param sigma `[numeric(1)=NULL]` \cr Standard deviation of the Gaussian
#' function. It related to the Gaussian decay rate \eqn{\lambda} as
#' `lambda = 1/(2*sigma^2)`. Only considered to compute the ZoI
#' for the `gaussian_decay` function when the ZoI radius parameter is null
#' (`radius = NULL`).
#'
#' @param round_vals `[numeric(1)=NULL]` \cr Number of digits for rounding the weights
#' in the output matrix. If `NULL` (default), weights are not rounded.
#' @param save_txt `[logical(1)=FALSE]` \cr Should the ZoI matrix be saved in an external
#' text file? If `FALSE` (default), the output matrix is just printed within the R session.
#'
#' @param save_format `[character(1)="GRASS_rmfilter"]{"GRASS_rmfilter", "raw"}` \cr
#' Format in which the function should be saved. Currently, either of the two options:
#' - GRASS GIS format for the module `r.mfilter`
#' (`save_format = "GRASS_rmfilter"`), see details [here](https://grass.osgeo.org/grass78/manuals/r.mfilter.html));
#' - raw matrix (`save_format = "raw"`), in which only the values of the matrix are printed.
#' @param save_folder `[character(1)=NULL]` \cr Path to the folder where the matrix file should be written.
#' If `NULL`, the current working directory is used.
#' @param save_file `[character(1)=NULL]` \cr Name of the output file, generally a ".txt" file.
#' If `NULL`, a standard filename is created, using the `type` and `radius`. E.g. "filter_bartlett2000.txt".
#' @param normalize `[logical(1)=FALSE]` \cr Whether the matrix should be normalized (sum of all cells is 1 if
#' `normalize = TRUE`) or kept as it is (default, `normalize = FALSE`).
#' @param divisor `[numeric(1)=1]` \cr By default, 1. This is the divisor of the neighborhood
#' matrix when used within `r.mfilter`. According the the module documentation, "The filter process produces a new
#' category value for each cell in the input raster map layer by multiplying the category values of the cells
#' in the n x n neighborhood around the center cell by the corresponding matrix value and adding them together.
#' If a divisor is specified, the sum is divided by this divisor." \cr
#' If the divisor is zero, "then the divisor is computed for each cell as the sum of the MATRIX values where
#' the corresponding input cell is non-null." In other words, the output map will be rescaled to the
#' interval $[0,1]$. If `normalize = TRUE`, the divisor is set to `n*n`.
#' @param parallel `[logical(1)=TRUE]` \cr Whether the computation should be paralelized or not (details in
#' the documentation of the [`r.mfilter`](https://grass.osgeo.org/grass78/manuals/r.mfilter.html) module).
#'
#' @return A matrix with the weight values. In the context of cumulative impact assessment, we call it a
#' zone of influence (ZoI) matrix used to compute the cumulative zone of influence. If `save_txt = TRUE`,
#' the matrix is saved in an output text file, e.g. to be used with external software.
#'
#' @example examples/filter_create_example.R
#'
#' @seealso See [oneimpact::zoi_functions()] for some ZoI function shapes and
#' [oneimpact::filter_save()] for options to save the ZoI matrix as a text file. \cr
#' See also [smoothie::kernel2dmeitsjer()], [terra::focalMat()], and
#' [raster::focalWeight()] for other functions to create filters or weight matrices. \cr
#' See
#' [r.mfilter](https://grass.osgeo.org/grass80/manuals/r.mfilter.html),
#' [r.resamp.filter](https://grass.osgeo.org/grass80/manuals/r.resamp.filter.html), and
#' [r.neighbors](https://grass.osgeo.org/grass80/manuals/r.neighbors.html) for
#' GRASS GIS uses of filters in neighborhood analysis.
#'
#' @export
filter_create <- function(r = 100,
                          radius = NULL,
                          type = c("exp_decay", "bartlett", "circle", "threshold_decay",
                                   "gaussian_decay", "Gauss", "rectangle")[1],
                          zoi_limit = 0.05,
                          half_life = NULL,
                          zoi_hl_ratio = NULL,
                          sigma = NULL,
                          min_intensity = 0.01,
                          max_dist = 5000,
                          normalize = FALSE,
                          divisor = 1,
                          round_vals = NULL,
                          save_txt = FALSE,
                          save_format = c("GRASS_rmfilter", "raw")[1],
                          save_folder = NULL,
                          save_file = NULL,
                          parallel = TRUE) {

  # check the input data class of r
  if(class(r) %in% c("RasterLayer", "RasterBrick", "RasterStack", "SpatRaster")) {
    res <- terra::res(r)[1]
  } else {
    if(is.numeric(r) & r > 0) {
      res <- r
    } else
      stop("'r' must be either an input raster map or a numeric value corresponding to the resolution of a raster.")

  }

  # apply function
  if(type == "exp_decay") {
    parms <- set_filt_exp_decay(radius = radius,
                                zoi_limit = zoi_limit,
                                res = res,
                                half_life = half_life,
                                zoi_hl_ratio = zoi_hl_ratio,
                                min_intensity = min_intensity,
                                max_dist = max_dist)
  }

  if(type %in% c("step", "threshold", "circle", "threshold_decay", "step_decay")) {
    parms <- set_filt_step(radius = radius, res = res)
  }

  if(type %in% c("bartlett", "batlett_decay", "tent_decay", "linear_decay")) {
    parms <- set_filt_bartlett(radius = radius, res = res)
  }

  if(type %in% c("rectangle", "box")) {
    parms <- set_filt_rectangle(radius = radius, res = res)
  }

  if(type %in% c("Gauss", "gauss", "gaussian", "gaussian_decay")) {
    parms <- set_filt_gassian_decay(radius = radius,
                                    zoi_limit = zoi_limit,
                                    res = res,
                                    sigma = sigma,
                                    min_intensity = min_intensity,
                                    max_dist = max_dist)
  }

  # get parameters
  radius <- parms$radius
  radius_pix <- parms$radius_pix
  size_pix <- parms$size_pix

  # create distance matrix
  # distance in pixels to the central cell of the matrix
  dist_mat <- sqrt((matrix(c(1:size_pix), nrow = size_pix, ncol = size_pix, byrow = F) - (radius_pix + 1))^2+
                     (matrix(c(1:size_pix), nrow = size_pix, ncol = size_pix, byrow = T) - (radius_pix + 1))^2)
  # image(dist_mat)
  # plot(terra::rast(dist_mat))

  # apply function
  if(type == "exp_decay") {
    dist_mat <- exp(-parms$lambda * dist_mat)
  }

  if(type %in% c("step", "threshold", "circle", "threshold_decay", "step_decay")) {
    dist_mat <- 1 * (dist_mat*res <= radius)
  }

  if(type %in% c("bartlett", "batlett_decay", "tent_decay", "linear_decay")) {
    dist_mat <- pmax((1 + parms$lambda * dist_mat), 0)
  }

  if(type %in% c("rectangle", "box")) {
    dist_mat[] <- 1
  }

  if(type %in% c("Gauss", "gauss", "gaussian", "gaussian_decay")) {
    dist_mat <- exp(-parms$lambda * dist_mat**2)
  }
  # image(dist_mat)
  # plot(terra::rast(dist_mat))

  # normalize
  if(normalize)
    # dist_mat <- dist_mat/sum(dist_mat[1+radius_pix,])
    dist_mat <- dist_mat/sum(dist_mat)

  # image(dist_mat)
  # plot(terra::rast(dist_mat))

  # round decimals
  if(!is.null(round_vals))
    if(round_vals >= 0) dist_mat <- round(dist_mat, round_vals)
  # image(dist_mat)
  # plot(terra::rast(dist_mat))

  if(save_txt) {
    # save matrix outside R for use within GRASS GIS
    oneimpact::filter_save(filt = dist_mat, radius = radius, type = type,
                           save_format = save_format, save_folder = save_folder,
                           save_file = save_file, parallel = parallel,
                           divisor = divisor, separator = " ")

  }

  dist_mat
}

set_filt_exp_decay <- function(radius = NULL,
                               zoi_limit = 0.05,
                               half_life = NULL,
                               res = 100,
                               zoi_hl_ratio = NULL,
                               min_intensity = 0.01,
                               max_dist = 50000){

  # define lambda depending on the input parameter
  if(!is.null(radius)) {

    # define radius in terms on number of pixels
    radius <- radius/res

    if(is.null(zoi_hl_ratio)) {
      lambda <- log(1/zoi_limit) / radius
    } else {
      half_life <- radius/zoi_hl_ratio
      lambda <- log(2)/half_life
    }

  } else {

    if(!is.null(half_life)) {
      # define radius or half life, depending on which is given as input
      half_life <- half_life/res
      lambda <- log(2)/half_life
    } else {
      stop("Either both 'radius' and 'zoi_limit' must be specified, or both 'half_life' and 'zoi_hl_ratio'.")
    }
  }

  # tmp <- exp(-lambda * c(0:round(half_life*6))/half_life)
  # define radius and size (diameter)
  tmp <- exp(-lambda * c(0:round(2*radius)))
  radius_pix <- min(which(tmp < min_intensity)[1], round(max_dist/res))
  size_pix <- 2*radius_pix + 1

  return(list(radius = radius, radius_pix = radius_pix, size_pix = size_pix, lambda = lambda))
}

set_filt_step <- function(radius, res){

  # define radius and size (diameter)
  radius_pix <- ceiling(radius/res)
  size_pix <- 2*radius_pix + 1

  return(list(radius = radius, radius_pix = radius_pix, size_pix = size_pix, lambda = NULL))
}

set_filt_rectangle <- function(radius, res){

  # define radius and size (diameter)
  radius_pix <- floor(radius/res)
  size_pix <- 2*radius_pix + 1

  return(list(radius = radius, radius_pix = radius_pix, size_pix = size_pix, lambda = NULL))
}

set_filt_bartlett <- function(radius, res){

  # define radius and size (diameter)
  radius_pix <- ceiling(radius/res)
  size_pix <- 2*radius_pix + 1
  # define beta (beta = -b/a or beta = -1/radius)
  lambda <- -1/(radius/res)

  return(list(radius = radius, radius_pix = radius_pix, size_pix = size_pix, lambda = lambda))
}

set_filt_gassian_decay <- function(radius = NULL,
                                   zoi_limit = 0.05,
                                   res = 100,
                                   sigma = NULL,
                                   min_intensity = 0.01,
                                   max_dist = 50000){


  if(!is.null(radius)) {
    # define radius in terms on number of pixels
    radius <- radius/res
    lambda = log(1/zoi_limit) / (radius**2)
  } else {
    if(!is.null(sigma)) {
      # define sigma in terms on number of pixels
      sigma <- sigma/res
      lambda = 1/(2*sigma**2)
    } else {
      stop("Either 'radius' or 'sigma' must be specified.")
    }

  }

  # tmp <- exp(-lambda * c(0:round(half_life*6))/half_life)
  # define radius and size (diameter)
  tmp <- exp(-lambda * c(0:round(2*radius))**2)
  radius_pix <- min(which(tmp < min_intensity)[1], round(max_dist/res))
  size_pix <- 2*radius_pix + 1

  return(list(radius = radius, radius_pix = radius_pix, size_pix = size_pix, lambda = lambda))
}
