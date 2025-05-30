#' Binarize continuous raster maps
#'
#' Through GRASS GIS, this function transforms continuous
#' or categorical maps (with more than 1 category)
#' into binary maps (with only two values), to represent, e.g.,
#' habitat-matrix maps in the context of landscape ecology. It can
#' also be used to produce binary maps form maps with only one category (and
#' all the rest as NULL/no-data). It requires an active connection between
#' the R session and a GRASS GIS location and mapset
#' (through the package [rgrass]), and that the input map is already loaded
#' within this GRASS GIS mapset.
#'
#' For a similar procedure within R, use raster algebra functions within the [raster] and [terra] packages.
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
#' @param overwrite `[logical(1)]` \cr Whether the output maps should be overwriten
#' (flag `overwrite = TRUE`).
#' @param input_as_region `[logical(1)=FALSE]` \cr Should the input map `x` be
#' used to redefine the working region in GRASS before raster binarization?
#' If `TRUE`, `x` is used to define the region with `g.region`. If `FALSE`,
#' the region previously defined in the GRASS GIS session is used for computation.
#' Default is `FALSE`.
#' @param verbose `[logical(1)=FALSE]` \cr Should messages of the computation steps
#' be printed in the prompt along the computation?
#'
#' @return A binarized map with only two values (or a set of binarized maps if `length(breaks)` > 1)
#' within the GRASS GIS mapset. In R, the output is a string with the name of this
#' map.
#'
#' @seealso See also the documentation of raster algebra with [terra] [here](https://rspatial.org/terra/pkg/4-algebra.html)
#' and with [raster] [here](https://rspatial.org/raster/pkg/4-algebra.html).
#'
#' @example examples/grass_binarize_example.R
#'
#' @export
grass_binarize <- function(x,
                           breaks = 0.5,
                           output = paste0(x, "_bin"),
                           null = NULL,
                           setnull = NULL,
                           bin_values = c(0, 1),
                           input_as_region = FALSE,
                           verbose = FALSE,
                           overwrite = FALSE, ...) {

  # Check function arguments ----
  if(is.numeric(breaks) == FALSE) stop("'breaks' must be a numeric vector.")

  # flags
  flags <- c()
  if(!verbose) flags <- c(flags, "quiet")
  if(overwrite) flags <- c(flags, "overwrite")

  # flags for g.region
  flags_region <- c("a")
  if(verbose) flags_region <- c(flags_region, "p")

  # output names
  out_name <- output
  if(length(breaks) > 1) out_name = paste0(out_name, "_", breaks)

  for(i in 1:length(breaks)) {

    # region
    if(input_as_region) {
      if(verbose) {
        msg  <- "Important: the GRASS GIS computational is being reset to the extent of the input map `x`."
        rgrass::execGRASS("g.message", message = msg)
      }

      rgrass::execGRASS("g.region", raster = x, flags = flags_region)
    }

    # nulls
    if(!(is.null(null) & is.null(setnull))) {

      # copy input
      inter_map = "inter_map" # intermediate map
      rgrass::execGRASS("g.copy",
                         parameters = list(raster = paste0(x, ",", inter_map)),
                         flags = flags)

      # set parameters
      parms <- list(map = inter_map)
      if(!is.null(null)) parms <- append(parms, list(null = null))
      if(!is.null(setnull)) parms <- append(parms, list(setnull = setnull))
      if(!verbose) qq <- "quiet" else qq <- c()
      # run r.null
      rgrass::execGRASS("r.null", parameters = parms, flags = qq)

      x <- inter_map
    }

    # binarize
    out <- out_name[i]

    expr = sprintf("if(A >= %f, %d, %d)", breaks[i], bin_values[2], bin_values[1])
    rgrass::execGRASS("r.mapcalc.simple",
                       expression = expr,
                       a = x,
                       output = out,
                       flags = flags)

  }

  # remove intermediate map
  if(!(is.null(null) & is.null(setnull))) {
    rgrass::execGRASS("g.remove", type = "raster", name = "inter_map",
                       flags = "f")
  }

  return(out_name)
}


