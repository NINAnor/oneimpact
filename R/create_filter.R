#' Create filters or kernel matrices for raster neighborhood analyses
#'
#' This function creates matrices of weights following different
#' functions to be used in neighborhood analyses for rasters. It is possible
#' to export these matrices as text files, for use with external softwares
#' such as the r.mfilter module within GRASS GIS.
#'
#' #### COMPARE WITH smoothie::kernel2dsmooth and smoothie::kernel2dmeitsjer, maybe wrap some options here
#' #### POLSSIBLY: IMPLEMENT IN THE SAME WAY AS FOCAL, WITH INPUT RASTER AS ARGUMENT, POSSIBLY
#' #### check: if the outer ring of the matrix is all zero, remove it
#'
#' @param r Either a numeric value corresponding to the resolution (pixel size) that each pixel in the filter matrix
#' should correspond to; or a raster object (`SpatRaster` from the `terra` package or `RasterLayer`, `RasterBrick`, or
#' `RasterStack` from the `raster` package) from which such resolution can be extracted.
#'
#' @param zoi `[numeric(1)=NULL]` \cr Zone of Influence (ZoI), in map units (preferentially meters).
#' The ZoI is the distance, scale, or buffer size around a feature up to which we consider there is
#' an effect or influence of an infrastructure or variable. \cr
#' For the circle neighborhood (equivalent to step or threshold neighborhoods), the `zoi` corresponds to the radius
#' (or thresold) of the circle, beyond which the filter is zero. \cr
#' For the rectangular neighborhood, the `zoi` corresponds to half the size of the square size, or
#' `square size = 2*zoi`. For a rectangular filter with different size of the sides, use [terra::focal()] (but
#' please note the interpretation of the parameters is different). \cr
#' For the Bartlett neighborhood, the `zoi` corresponds to the distance beyond which the filter is zero. \cr
#' For the exponential decay neighborhood, the `zoi` is used to define the half-life and the lambda
#' of the exponential decay function, based on the parameter `zoi_hl_ratio`,
#' that defines the ratio between the ZoI and the exponential half-life. Since the half-life is the value
#' where the exponential decay decreases by `0.5`, a ratio of, for instance, `zoi_hl_ratio = 4` (default)
#' would mean that the ZoI is defined as the value where the exponential decay decreases to `0.5^4 = 0.0625`.
#' In this case, if `zoi = 4000` m, this means that the ZoI is four times higher than the half-life, i.e.
#' `half_life = 1000` and `lambda = log(2)/half_life = 6.93e-4`. The definition of a zone of
#' influence does not imply a cutoff of the exponential decay function but is only used to define
#' its parameters, based on the defined `zoi_hl_ratio` parameter. The cutoff is given either by the
#' `min_intensity` or the `max_dist` parameters.
#' If `zoi = NULL`, the exponential decay is defined based on the `half_life` parameter. \cr
#' Gaussian neighborhood with parameterization in terms of `zoi` to be implemented.
#'
#' @param method `[character(1)="exp_decay"]{"exp_decay", "bartlett", "circle", "threshold", "step", "rectangle"}` \cr
#' Rectangle = boxcar in smoothie::kernel2dmeitsjer
#' Gaussian neighborhood with parameterization in terms of `zoi` to be implemented. So far it is possible to
#' create Gaussian filters using other functions such as [terra::focalMat()] with parameter
#' `type = "Gauss"` and [smoothie::kernel2dmeitsjer()] with parameter `type = "gauss"`.
#'
#' @param half_life `[numeric(1)=NULL]` \cr Half life parameter of the exponential decay function, in meters. If NULL,
#' the half life is define in terms of the ZoI and the `zoi_hl_ratio` parameter, which defines the ratio
#' between the ZoI and the half life. By default, we set this ratio as `zoi/half_life = 4`.
#' The exponent of the exponential decay distance function is defined as `lambda = log(2)/half_life`.
#' @param zoi_hl_ratio `[numeric(1)=6]` \cr Ratio between the ZoI and the half life of the exponential decay
#' distance function. It is used to define the ZoI for the exponential decay function. For instance, if
#' `half_life = 1000` and `zoi_hl_ratio = 4`, the ZoI will be 4000 m (when the exponential decay decrease to
#' `0.5**4 = 0.0625`.
#' @param min_intensity `[numeric(1)=0.01]` \cr Minimum intensity of the exponential decay function to
#' define the size (radius) of the window that define the filter.
#' @param max_dist `[numeric(1)=50000]` \cr Maximum size (in meters) to
#' define the size (radius) of the window that define the filter.
#'
#' @param divisor ...
#' @param save_format `[character(1)="GRASS_r.mfilter"]{"GRASS_r.mfilter", "raw"}` \cr
#' Format in which the function should be saved. Currently, only GRASS GIS format for the module `r.mfilter`
#' (`save_format = "GRASS_r.mfilter"`, according to the required format for `r.mfilter` module, details
#' [here](https://grass.osgeo.org/grass78/manuals/r.mfilter.html)) or raw (`save_format = "raw"`),
#' in which only the values of the matrix are printed.
#' @param save_folder `[character(1)=NULL]` \cr Path to the folder where the matrix file should be written.
#' If `NULL`, the current directory is used.
#' @param save_file `[character(1)=NULL]` \cr Name of the output file, generally a ".txt" file.
#' If `NULL`, a standard filename is created, using the the `method` and `zoi`. E.g. "filter_bartlett2000.txt".
#' @param normalize `[logical(1)=FALSE]` \cr Whether the matrix should be normalized (sum of all cell is 1, if
#' `normalize = TRUE`) or kept as it is (default, `normalize = FALSE`).
#' @param divisor `[numeric(1)=1]` \cr By default, 1. This is the divisor of the neighborhood
#' matrix, when used within `r.mfilter`. According the the module documentation, "The filter process produces a new ¨
#' category value for each cell in the input raster map layer by multiplying the category values of the cells
#' in the n x n neighborhood around the center cell by the corresponding matrix value and adding them together.
#' If a divisor is specified, the sum is divided by this divisor." \cr
#' If the divisor is zero, "then the divisor is computed for each cell as the sum of the MATRIX values where
#' the corresponding input cell is non-null." In other words, the output map will be rescaled to the
#' interval [0,1]. If `normalize = TRUE`, the divisor is set to `n*n`.
#' @param parallel `[logical(1)=TRUE]` \cr Whether the computation should be paralelized or not (details in
#' the documentation of the [`r.mfilter`]((https://grass.osgeo.org/grass78/manuals/r.mfilter.html)) module).
#'
#' @return A matrix with the weight values.
#'
#' @example examples/create_filter_example.R
#'
#' @seealso See also [smoothie::kernel2dmeitsjer()], [terra::focalMat()], and
#' [raster::focalWeight()] for other functions to create filters or weight matrices.
#'
#' @export
create_filter <- function(r = 100,
                          zoi = NULL,
                          method = c("exp_decay", "bartlett", "circle", "threshold", "step", "Gauss", "rectangle")[1],
                          half_life = NULL,
                          zoi_hl_ratio = 4,
                          min_intensity = 0.01,
                          max_dist = 5000,
                          normalize = FALSE,
                          divisor = 1,
                          round_vals = NULL,
                          save_txt = FALSE,
                          save_format = c("GRASS_r.mfilter", "raw")[1],
                          save_folder = NULL,
                          save_file = NULL,
                          output = c("CumInf", "Densiy")[1],
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
  if(method == "exp_decay") {
    parms <- set_filt_exp_decay(zoi = zoi, half_life = half_life,
                                res = res, zoi_hl_ratio = zoi_hl_ratio,
                                min_intensity = min_intensity,
                                max_dist = max_dist)
  }

  if(method %in% c("step", "threshold", "circle")) {
    parms <- set_filt_step(zoi = zoi, res = res)
  }

  if(method == "bartlett") {
    parms <- set_filt_bartlett(zoi = zoi, res = res)
  }

  if(method == "rectangle") {
    parms <- set_filt_rectangle(zoi = zoi, res = res)
  }

  # get parameters
  zoi <- parms$zoi
  radius_pix <- parms$radius_pix
  size_pix <- parms$size_pix

  # create distance matrix
  # distance in pixels to the central cell of the matrix
  dist_mat <- sqrt((matrix(c(1:size_pix), nrow = size_pix, ncol = size_pix, byrow = F) - (radius_pix + 1))^2+
                     (matrix(c(1:size_pix), nrow = size_pix, ncol = size_pix, byrow = T) - (radius_pix + 1))^2)
  # image(dist_mat)
  # plot(terra::rast(dist_mat))

  # apply function
  if(method == "exp_decay") {
    dist_mat <- exp(-parms$lambda * dist_mat)
  }

  if(method %in% c("step", "threshold", "circle")) {
    dist_mat <- 1 * (dist_mat*res <= zoi)
  }

  if(method == "bartlett") {
    dist_mat <- pmax((1 + parms$lambda * dist_mat), 0)
  }

  if(method == "rectangle") {
    dist_mat[] <- 1
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
     if(round_vals > 0) dist_mat <- round(dist_mat, round_vals)
  # image(dist_mat)
  # plot(terra::rast(dist_mat))

  if(save_txt) {
    # save matrix outside R for use within GRASS GIS
    save_mfilter(filt = dist_mat, zoi = zoi, method = method,
                 save_format = save_format, save_folder = save_folder,
                 save_file = save_file, parallel = parallel,
                 divisor = divisor, separator = " ")

  }

  dist_mat
}

set_filt_exp_decay <- function(zoi = NULL,
                               half_life = NULL,
                               res = 100,
                               zoi_hl_ratio = 4,
                               min_intensity = 0.01,
                               max_dist = 50000){

  # define zoi or half life, depending on which is given as input
  if(is.null(zoi))
    zoi <- half_life * zoi_hl_ratio
  else
    half_life <- zoi/zoi_hl_ratio

  # define half life in terms on number of pixels
  half_life <- half_life/res
  # define lambda
  lambda <- log(2)/half_life
  # tmp <- exp(-lambda * c(0:round(half_life*6))/half_life)
  # define radius and size (diameter)
  tmp <- exp(-lambda * c(0:round(2*zoi)))
  radius_pix <- min(which(tmp < min_intensity)[1], round(max_dist/res))
  size_pix <- 2*radius_pix + 1

  return(list(zoi = zoi, radius_pix = radius_pix, size_pix = size_pix, lambda = lambda))
}

set_filt_step <- function(zoi, res){

  # define radius and size (diameter)
  radius_pix <- ceiling(zoi/res)
  size_pix <- 2*radius_pix + 1

  return(list(zoi = zoi, radius_pix = radius_pix, size_pix = size_pix, lambda = NULL))
}

set_filt_rectangle <- function(zoi, res){

  # define radius and size (diameter)
  radius_pix <- floor(zoi/res)
  size_pix <- 2*radius_pix + 1

  return(list(zoi = zoi, radius_pix = radius_pix, size_pix = size_pix, lambda = NULL))
}

set_filt_bartlett <- function(zoi, res){

  # define radius and size (diameter)
  radius_pix <- ceiling(zoi/res)
  size_pix <- 2*radius_pix + 1
  # define beta (beta = -b/a or beta = -1/zoi)
  lambda <- -1/(zoi/res)

  return(list(zoi = zoi, radius_pix = radius_pix, size_pix = size_pix, lambda = lambda))
}
