#' Binarize continuous raster maps
#'
#' Through GRASS GIS, this function transforms continuous
#' or categorical maps (with more than 1 category)
#' into binary maps (with only two values), to represent, e.g.,
#' habitat-matrix maps in the context of landscape ecology. It can
#' also be used to produce binary maps form maps with only one category (and
#' all the rest as NULL/no-data). It requires an active connection between
#' the R session and a GRASS GIS location and mapset
#' (through the package [rgrass7]), and that the input map is already loaded
#' within this GRASS GIS mapset.
#'
#' For a similar procedure within R, use the function [landscapetools::util_binarize()]
#' or the raster algebra functions within the [raster] and [terra] packages.
#'
#' @param x `[character(1)]` \cr Name of the input raster map, within a GRASS GIS location and mapset.
#'
#' @param breaks `[numeric]` \cr Breaks or threshold to define the binary values in the output
#' binary map. All pixels with `value >= breaks` are considered as 1 (or the upper value defined in
#' `bin_values`), and all the rest are considered as 0 (or the lower value defined in `bin_values`).
#' `breaks` might be either a single numeric value or a vector of numeric values, in which case
#' multiple binary maps are created (with different break thresholds).
#'
#' @param output `[character(1)]` \cr Name of the output map, or prefix of the output map if
#' `length(breaks) > 1`. In the latter case, the names are completed with the break value.
#' The defult is to use the same name as the input map, plus "_bin" in the end.
#'
#' @param null `[numeric(1)=NULL]` \cr If `NULL` (default), all NULL/no-data pixels in from `x`
#' are kept as they are in the output map. Otherwise, a numeric value that all NULL pixels should
#' assume in the output map. It uses the module
#' [r.null](https://grass.osgeo.org/grass78/manuals/r.null.html)).
#'
#' @param setnull `[]` \cr If `NULL` (default), no changes are made. Otherwise, a set of numeric
#' values that should be transformed into NULL/NA data (using the module
#' [r.null](https://grass.osgeo.org/grass78/manuals/r.null.html)).
#' @param bin_values `[numeric(2)=c(0,1)]` \cr Values c(lower, upper) that the output map pixels should
#' have if their values are either "lower" or "equal or higher" `breaks`. By default, c(0, 1).
#' @param quiet,overwrite `[logical(1)]` \cr Whether the procedures is GRASS GIS should be run
#' quetly (flag `quiet = TRUE`) and whether the output maps should be overwriten (flag `overwrite = TRUE`).
#'
#' @return A binarized map with only two values (or a set of binarized maps if `length(breaks)` > 1).
#'
#' @seealso See also [landscapetools::util_binarize()], [landscapetools::util_classify()],
#' and a documentation of raster algebra with [terra] [here](https://rspatial.org/terra/pkg/4-algebra.html) and with
#' [raster] [here](https://rspatial.org/raster/pkg/4-algebra.html).
#'
#' @example examples/util_binarize_GRASS_example.R
#'
#' @export
util_binarize_GRASS <- function(x,
                                breaks = 0.5,
                                output = paste0(x, "_bin"),
                                null = NULL,
                                setnull = NULL,
                                bin_values = c(0, 1),
                                quiet = TRUE,
                                overwrite = FALSE, ...) {

  # Check function arguments ----
  if(is.numeric(breaks) == FALSE) stop("'breaks' must be a numeric vector.")

  # flags
  flags <- c()
  if(quiet) flags <- c(flags, "quiet")
  if(overwrite) flags <- c(flags, "overwrite")

  # flags for g.region
  flags_region <- c("a")
  if(!quiet) flags_region <- c(flags_region, "p")

  # output names
  out_name <- output
  if(length(breaks) > 1) out_name = paste0(out_name, "_", breaks)

  for(i in 1:length(breaks)) {

    # region
    rgrass7::execGRASS("g.region", raster = x, flags = flags_region)

    # binarize
    out <- out_name[i]

    expression = sprintf("if(A >= %f, %d, %d)", breaks[i], bin_values[2], bin_values[1])
    rgrass7::execGRASS("r.mapcalc.simple",
                       expression = expression,
                       a = x,
                       output = out,
                       flags = flags)

    if(!(is.null(null) & is.null(setnull))) {
      # set parameters
      parms <- list(map = out)
      if(!is.null(null)) parms <- append(parms, list(null = null))
      if(!is.null(setnull)) parms <- append(parms, list(setnull = setnull))
      if(quiet) qq <- "quiet" else qq <- c()
      # run r.null
      rgrass7::execGRASS("r.null", parameters = parms, flags = qq)
    }
  }

  return(out_name)
}


