#' Rasterizes a vector counting the number of features in each pixel
#'
#' Add within and other arguments for grass functions as options/parameters
#'
#' @param x Input vector map.
#' @param output Output map name.
#' @param column `[chracter(1)=NULL` \cr Default is `NULL`. If not `NULL`, the name of a column
#' in the input vector `x` that corresponds to the column to be summed to count the number
#' of features in each pixel of the output raster map. If `NULL`, this column is created in
#' a temporary vector, with all values equal 1.
#' @param set_region `[logical(1)=TRUE]` \cr Default is TRUE. Whether the GRASS GIS
#' computational region should be set within the function (to the extent of `x`) or not.
#' If `FALSE`, the current computational region is used.
#' @param align `[character(1)=NULL]` \cr Name of a raster map with which to align the
#' computational region to produce the output map.
#'
#' @example examples/util_v2rast_count_GRASS_example.R
#'
#' @export
util_v2rast_count_GRASS <- function(x,
                                    output = paste0(x, "_count"),
                                    column = NULL,
                                    set_region = TRUE,
                                    align = NULL,
                                    remove_intermediate = TRUE,
                                    quiet = TRUE,
                                    verbose = FALSE,
                                    overwrite = FALSE, ...) {

  #--- options ---
  # flags
  flags <- c()
  flags_text <- ""
  if(quiet) flags <- c(flags, "quiet")
  if(overwrite) flags <- c(flags, "overwrite")

  # flags for g.region
  flags_region <- c("a")
  if(verbose) flags_region <- c(flags_region, "p")

  # intermediate maps to remove
  if(remove_intermediate) to_remove <- c()

  #--- region ---
  # region
  if(set_region)
    if(is.null(align)) {
      rgrass7::execGRASS("g.region", vector = x, flags = flags_region)
    } else
      rgrass7::execGRASS("g.region", vector = x, align = align, flags = flags_region)

  #--- procedures ---
  # set temporary region
  temp_region <- "temp_region_v2rast_count"
  if(remove_intermediate) to_remove <- c(to_remove, temp_region)
  rgrass7::execGRASS("v.in.region", output = temp_region, flags = flags)

  # set temporary vector
  if(verbose) rgrass7::execGRASS("g.message", message = "Copying input vector...")
  temp_vect <- "temp_vect_v2rast_count"
  if(remove_intermediate) to_remove <- c(to_remove, temp_vect)
  flags_extract <- c("t", flags)
  rgrass7::execGRASS("v.select", ainput = x, binput = temp_region,
                     operator = "within", output = temp_vect, flags = flags_extract)

  # add new column for counting
  if(verbose) rgrass7::execGRASS("g.message", message = "Adding new column for counting...")
  if(is.null(column)) {
    column_count <- "val_c"
  } else {
    column_count <- column
  }
  flags_table <- c()
  if(quiet) flags_table <- c(flags_table, "quiet")
  rgrass7::execGRASS("v.db.addtable", map = temp_vect, columns = sprintf("%s integer", column_count),
                     flags = flags_table)
  rgrass7::execGRASS("v.db.update", map = temp_vect, column = column_count, value = "1",
                     flags = flags_table)

  # write temporary ascii file
  if(verbose) rgrass7::execGRASS("g.message", message = "Writing output ascii and rasterizing...")
  f <- tempfile(fileext = ".txt")
  rgrass7::execGRASS("v.out.ascii", input = temp_vect, output = f, columns = column_count,
                     flags = flags)
  # read it as raster
  rgrass7::execGRASS("r.in.xyz", input = f, z = 4, output = output, method = "sum",
                     flags = flags)

  # remove intermediate maps
  remove_flags = ifelse(quiet, c("f", "quiet"), "f")
  if(remove_intermediate) rgrass7::execGRASS("g.remove", type = "vect", name = to_remove,
                                             flags = remove_flags)

  # return output map names
  output
}
