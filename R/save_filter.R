#' Save kernel/filter matrix to use in r.mfilter within GRASS GIS
#'
#' This function saves a matrix with weights (filter or kernel matrix) in an external text file.
#' It can save either the raw matrix
#' or use the standards for running `r.mfilter` within GRASS GIS
#' (with specific header and details).
#'
#' @details If used in the `r.mfilter` GRASS GIS module, "The filter process produces a new category value
#' for each cell in the input raster map layer by multiplying the category values of the cells in the n x n
#' neighborhood around the center cell by the corresponding matrix value and adding them together.
#' If a divisor is specified, the sum is divided by this divisor." See details
#' [here](https://grass.osgeo.org/grass78/manuals/r.mfilter.html).
#'
#' @param filt `[matrix]` \cr Filter or weight matrix, such as one created by [oneimpact::create_filter] or [terra::focalMat].
#' @param zoi_radius `[numeric(1)]` \cr Radius of the Zone of Influence (ZoI) of the matrix, in meters.
#' @param type `[character(1)]` Function for the kernel or filter matrix (see `type` parameter for [oneimpact::create_filter]).
#' @param save_format `[character(1)="GRASS_rmfilter"]{"GRASS_rmfilter", "raw"}` \cr
#' Format in which the function should be saved. Currently, either of the two options:
#' - GRASS GIS format for the module `r.mfilter`
#' (`save_format = "GRASS_rmfilter"`), see details [here](https://grass.osgeo.org/grass78/manuals/r.mfilter.html));
#' - raw matrix (`save_format = "raw"`), in which only the values of the matrix are printed.
#' @param save_folder `[character(1)=NULL]` \cr Path to the folder where the matrix file should be written.
#' If `NULL`, the current working directory is used.
#' @param save_file `[character(1)=NULL]` \cr Name of the output file, generally a ".txt" file.
#' If `NULL`, a standard filename is created, using the `type` and `zoi_radius`. E.g. "filter_bartlett2000.txt".
#' @param normalize `[logical(1)=FALSE]` \cr Whether the matrix should be normalized (sum of all cells is 1 if
#' `normalize = TRUE`) or kept as it is (default, `normalize = FALSE`).
#' @param divisor `[numeric(1)=1]` \cr By default, 1. This is the divisor of the neighborhood
#' matrix when used within `r.mfilter`. According the the module documentation, "The filter process produces a new
#' category value for each cell in the input raster map layer by multiplying the category values of the cells
#' in the n x n neighborhood around the center cell by the corresponding matrix value and adding them together.
#' If a divisor is specified, the sum is divided by this divisor." \cr
#' If the divisor is zero, "then the divisor is computed for each cell as the sum of the MATRIX values where
#' the corresponding input cell is non-null." In other words, the output map will be rescaled to the
#' interval [0,1]. If `normalize = TRUE`, the divisor is set to `n*n`.
#' @param parallel `[logical(1)=TRUE]` \cr Whether the computation should be paralelized or not (details in
#' the documentation of the [`r.mfilter`](https://grass.osgeo.org/grass78/manuals/r.mfilter.html) module).
#' @param separator `[character(1)=" "]` \cr Separator between values of the matrix, within each line. Default is
#' a space.
#'
#' @return None. The funcion only saves the input matrix as an external file.
#'
#' @seealso [r.mfilter](https://grass.osgeo.org/grass78/manuals/r.mfilter.html)
#'
#' @examples
#' my_filter <- create_filter(r = 100, type = "bartlett", zoi_radius = 1000, round = 4)
#' save_filter(my_filter, zoi_radius = 1000, type = "bartlett", save_format = "GRASS_rmfilter")
#'
#' @export
save_filter <- function(
  filt,
  zoi_radius,
  type,
  save_format = c("GRASS_rmfilter", "raw")[1],
  save_folder = NULL,
  save_file = NULL,
  divisor = 1,
  normalize = FALSE,
  parallel = TRUE,
  separator = " ") {

  # define the divisor
  DIV <- divisor
  # if normalize = TRUE, replace DIV
  if(normalize) DIV <- size_pix**2

  # should the computation be parallelized?
  if(parallel) TYP <- "P" else TYP <- "S"

  # output file name
  if(is.null(save_file)) {
    save_file <- paste0("filter_", type, zoi_radius, ".txt")
  }
  if(is.null(save_folder)) {
    file_out <- save_file
  } else {
    file_out <- paste0(save_folder, "/", save_file)
  }

  # open file
  con <- file(file_out, "w")

  # r.mfilter preamble
  if(save_format == "GRASS_rmfilter") {
    writeLines(paste0("TITLE filter ", type, " ", zoi_radius, "m"), con = con, sep = "\n", useBytes = FALSE)
    writeLines(paste0("MATRIX ", dim(filt)[2]), con = con, sep = "\n", useBytes = FALSE)
  }

  # write matrix
  for (i in c(1:nrow(filt))) {
    writeLines(paste0(filt[i,], collapse = separator), con = con, sep = "\n", useBytes = FALSE)
  }

  # r.mfilter post-info
  if(save_format == "GRASS_rmfilter") {
    writeLines(paste0("DIVISOR ", DIV), con = con, sep = "\n", useBytes = FALSE)
    writeLines(paste0("TYPE ", TYP), con = con, sep = "\n", useBytes = FALSE)
  }

  # close file
  close(con)
}
