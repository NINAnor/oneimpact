#' Rasterizes a vector counting the number of features in each pixel
#'
#' This function rasterizes a vector file in GRASS GIS, counting the number of
#' vector features within each pixel. The function uses the resolution and extent
#' already set as the GRASS GIS mapset's computational region. If `input_as_region`
#' is set to `TRUE`, the extent of the vector map `x` is used reset as the region.
#' The resolution, though, continues the same set up earlier through `g.region`.
#'
#' @param x Input vector map.
#' @param output Output map name.
#' @param column `[chracter(1)=NULL` \cr Default is `NULL`. If not `NULL`, the name of a column
#' in the input vector `x` that corresponds to the column to be summed to count the number
#' of features in each pixel of the output raster map. If `NULL`, this column is created in
#' a temporary vector, with all values equal 1.
#' @param input_as_region `[logical(1)=FALSE]` \cr Default is FALSE. Whether the GRASS GIS
#' computational region should be set within the function (to the extent of `x`) or not.
#' If `FALSE`, the current computational region is used.
#' @param align `[character(1)=NULL]` \cr Name of a raster map with which to align the
#' computational region to produce the output map.
#' @param verbose `[logical(1)=FALSE]` \cr Should messages of the computation steps
#' be printed in the prompt along the computation?
#' @param overwrite `[logical(1)]` \cr Whether the output maps should be overwriten
#' (flag `overwrite = TRUE`).
#'
#' @return A raster map with the count of features within each pixel. The map is written
#' within the GRASS GIS mapset. In R, the output is a string with the name of this
#' map.
#'
#' @example examples/grass_v2rast_count_example.R
#'
#' @export
grass_v2rast_count <- function(x,
                               output = paste0(x, "_count"),
                               column = NULL,
                               input_as_region = FALSE,
                               align = NULL,
                               remove_intermediate = TRUE,
                               verbose = FALSE,
                               overwrite = FALSE, ...) {

  #--- options ---
  # flags
  flags <- c()
  flags_text <- ""
  if(!verbose) flags <- c(flags, "quiet")
  if(overwrite) flags <- c(flags, "overwrite")

  # flags for g.region
  flags_region <- c("a")
  if(verbose) flags_region <- c(flags_region, "p")

  # intermediate maps to remove
  if(remove_intermediate) to_remove <- c()

  #--- region ---
  # region
  if(input_as_region) {
    if(verbose) {
      msg  <- "Important: the GRASS GIS computational is being reset to the extent of the input map `x`."
      rgrass::execGRASS("g.message", message = msg)
    }

    if(is.null(align)) {
      rgrass::execGRASS("g.region", vector = x, flags = flags_region)
    } else
      rgrass::execGRASS("g.region", vector = x, align = align, flags = flags_region)
  }

  #--- procedures ---
  # If no column is given, the vector is copied only for the computational region and
  # a column with value=1 is created
  if(is.null(column)) {

    # set temporary region
    temp_region <- "temp_region_v2rast_count"
    if(remove_intermediate) to_remove <- c(to_remove, temp_region)
    rgrass::execGRASS("v.in.region", output = temp_region, flags = flags)

    # set temporary vector
    if(verbose) rgrass::execGRASS("g.message", message = "Copying input vector...")
    temp_vect <- "temp_vect_v2rast_count"
    if(remove_intermediate) to_remove <- c(to_remove, temp_vect)
    flags_extract <- c("t", flags)
    rgrass::execGRASS("v.select", ainput = x, binput = temp_region,
                       operator = "within", output = temp_vect, flags = flags_extract)

    # add new column for counting
    if(verbose) rgrass::execGRASS("g.message", message = "Adding new column for counting...")
    column_count <- "val_c"
    flags_table <- c()
    if(!verbose) flags_table <- c(flags_table, "quiet")
    rgrass::execGRASS("v.db.addtable", map = temp_vect, columns = sprintf("%s integer", column_count),
                       flags = flags_table)
    rgrass::execGRASS("v.db.update", map = temp_vect, column = column_count, value = "1",
                       flags = flags_table)
  } else {
    # If the column is given, the input vector and column are used
    temp_vect <- x
    column_count <- column
  }

  # write temporary ascii file
  if(verbose) rgrass::execGRASS("g.message", message = "Writing output ascii and rasterizing...")
  f <- tempfile(fileext = ".txt")
  rgrass::execGRASS("v.out.ascii", input = temp_vect, output = f, columns = column_count,
                     flags = flags)
  # read it as raster
  rgrass::execGRASS("r.in.xyz", input = f, z = 4, output = output, method = "sum",
                     flags = flags)

  # remove intermediate maps
  remove_flags = ifelse(!verbose, c("f", "quiet"), "f")
  if(remove_intermediate & length(to_remove) > 0)
    rgrass::execGRASS("g.remove", type = "vect", name = to_remove,
                       flags = remove_flags)

  # return output map names
  output
}
