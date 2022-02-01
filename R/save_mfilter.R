#' Save kernel/filter matrix to use in r.mfilter within GRASS GIS
#' 
#' @param filt `[matrix]` \cr Matrix, such as one created by [oneimpact::create_filter()] or [terra::focalMat()].
#' @param zoi `[numeric(1)]` \cr Zone of Influence of the matrix, in meters.
#' @param method `[character(1)]` Function for the kernel or filter matrix (see `method` parameter for [create_filter()]).
#' @param save_format `[character(1)=GRASS]{"GRASS", "raw"}` \cr Format in which the function should be saved. Currently, only GRASS GIS format 
#' (`save_format = "GRASS"`, according to the required format for `r.mfilter` module, details 
#' [here](https://grass.osgeo.org/grass78/manuals/r.mfilter.html)) or raw (`save_format = "raw"`),
#' in which only the values of the matrix are printed.
#' @param save_folder `[character(1)=NULL]` \cr Path to the folder where the matrix file should be written.
#' If `NULL`, the current directory is used.
#' @param save_file `[character(1)=NULL]` \cr Name of the output file, generally a ".txt" file. 
#' If `NULL`, a standard filename is created, using the the `method` and `zoi`. E.g. "filter_bartlett2000.txt".
#' @param normalize `[logical(1)=FALSE]` \cr Whether the matrix should be normalized (sum of all cell is 1, if
#' `normalize = TRUE`) or kept as it is (default, `normalize = FALSE`).
#' @param parallel `[logical(1)=TRUE]` \cr Whether the computation should be paralelized or not (details in 
#' the documentation of the [`r.mfilter`]((https://grass.osgeo.org/grass78/manuals/r.mfilter.html)) module).
#' @param separator `[character(1)=" "]` \cr Separator between values of the matrix, within each line. Default is 
#' a space.
#' 
#' @return None. The funcion only saves the input matrix as an external file.
#' 
#' @examples 
#' my_filter <- create_filter(r = 100, method = "bartlett", zoi = 1000, round = 4)
#' save_mfilter(my_filter, zoi = 1000, method = "bartlett", save_format = "GRASS")
#' 
#' @export
save_mfilter <- function(
  filt,
  zoi, 
  method, 
  save_format = c("GRASS", "raw")[1],
  save_folder = NULL,
  save_file = NULL,
  normalize = FALSE,
  parallel = TRUE,
  separator = " ") {
  
  if(!normalize) DIV <- 0 else DIV <- size_pix**2
  if(parallel) TYP <- "P" else TYP <- "S"
  if(save_format == "GRASS") {
    if(is.null(save_file)) {
      save_file <- paste0("filter_", method, zoi, ".txt")
    }
    if(is.null(save_folder)) { 
      file_out <- save_file
    } else {
      file_out <- paste0(save_folder, "/", save_file)
    } 
    
    con <- file(file_out, "w")
    writeLines(paste0("TITLE filter ", method, " ", zoi, "m"), con = con, sep = "\n", useBytes = FALSE)
    writeLines(paste0("MATRIX ", dim(filt)[2]), con = con, sep = "\n", useBytes = FALSE)
    for (i in c(1:nrow(filt))) {
      writeLines(paste0(filt[i,], collapse = separator), con = con, sep = "\n", useBytes = FALSE)
    }
    writeLines(paste0("DIVISOR ", DIV), con = con, sep = "\n", useBytes = FALSE)
    writeLines(paste0("TYPE ", TYP), con = con, sep = "\n", useBytes = FALSE)
    close(con)
  }
}