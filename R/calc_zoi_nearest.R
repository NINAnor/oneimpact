#' Calculate the zone of influence from the nearest feature
#'
#' @description
#' This function takes in a raster with locations or counts of
#' infrastructure and calculates a raster (or set of rasters, in case there is
#' more the one value for `radius`) representing the zone of influence (ZoI)
#' from the neareast feature of that type of infrastructure. Zones of influence
#' are defined by functions that decay with the distance from each
#' infrastructure and their rate of decay is controlled by the ZoI radius
#' (`radius`), which defines how far the influence of an infrastructure
#' feature goes. By default, the Gaussian decay ZoI is calculated, but other
#' decay shapes might be used (see [oneimpact::zoi_functions()] for examples).
#' The function might also return the distance to the nearest feature
#' or a transformation from it (e.g. log- and sqrt-distance from the nearest
#' feature).
#'
#' The procedure might be computed in both R and GRASS GIS. In GRASS, it
#' requires an active connection between the R session and a GRASS GIS
#' location and mapset (through the package [rgrass]), and that the input
#' maps are already loaded within this GRASS GIS mapset.
#' If the calculations are done in R, the input is a (set of) raster map(s)
#' and the function returns another (set of) raster map(s). If the calculations
#' are done within GRASS GIS, the input is the name of a raster map already
#' loaded in a GRASS GIS location and mapset, and the function returns
#' only the name of the output map. This map is stored in the GRASS GIS
#' location/mapset, and might be retrieved to R through the
#' [rgrass::read_RAST()] function or exported outside GRASS using the
#' `r.out.gdal` module, for instance.
#'
#' @details
#' In practice, the function `calc_zoi_nearest()` first calculates the distance from each
#' pixel to the nearest feature and then transforms it according to the ZoI
#' functions. In R, `calc_zoi_nearest()` makes use of the [terra::distance()]
#' function and the following procedures are made through raster algebra.
#' In GRASS, the module
#' [`r.grow.distance`](https://grass.osgeo.org/grass78/manuals/r.grow.distance.html)
#' is used to calculate the Euclidean distance from
#' the nearest feature and
#' [`r.mapcalc.simple`](https://grass.osgeo.org/grass78/manuals/r.mapcalc.simple.html)
#' to transform the distance into the different ZoI of the nearest feature.
#'
#' The input raster `x` should have positive values in the pixels where
#' infrastructure are located and NA/no-data in all other places.
#' The input raster is supposed to
#' represent the location of point, line, or polygon infrastructure
#' (e.g. houses, roads, mining areas), but any landscape variable whose
#' representation might be one of those would fit here
#' (e.g. areas of forest or any other habitat type or land cover).
#' We recommend that the input raster has a metric projection, so that distances
#' and zones of influence are based on distance to infrastructure measured in meters.
#'
#' @param x `[RasterLayer,SpatRaster]` \cr Raster representing locations of features,
#' preferentially with positive value where the features
#' are located and NA elsewhere. Alternatively, `x` might be a binary (dummy)
#' spatial variable representing the presence of linear or area features, with
#' NA/no-data as background.
#' `x` can be a `RasterLayer` from [raster] package or a [SpatRaster] from
#' [terra] package. If `where = "GRASS"`, `x` must be a string corresponding
#' to the name of the input map within a GRASS GIS location and mapset.
#'
#' The input raster `x` should have positive values in the pixels where infrastructure
#' are located and NA/no-data in all other places. In R it is also possible to have zeros
#' as the background and set `zeroAsNA = TRUE` for the computation of the ZoI of the
#' nearest feature.
#' In GRASS, maps without NA as background might be prepared as input for `calc_zoi_nearest()`
#' through [raster algebra](https://grass.osgeo.org/grass78/manuals/r.mapcalc.html)
#' and e.g. through the use of the module
#' [`r.null`](https://grass.osgeo.org/grass80/manuals/r.null.html).
#'
#' @param radius `[numeric(1)]` \cr Radius of the zone of influence (ZoI),
#' the distance at which the ZoI vanishes or goes below a given minimum limit value
#' `zoi_limit`. See [oneimpact::zoi_functions()] for details.
#' It can be a single value or a vector of values, in which case
#' several ZoI layers (one for each radius) are created. This parameter is
#' ignored if `type = "euclidean"`, `type = "log"`, or `type = "sqrt"`.
#'
#' @param type `[character(1)="Gauss"]{"Gauss", "exp_decay", "bartlett",
#' "threshold", "step", "euclidean", "log","sqrt"}` \cr
#' Shape of the zone of influence: \cr
#' \itemize{
#'   \item If `Gauss` or `half_norm`, the ZoI follows a half-normal shape: \cr
#'   `intercept * exp(-lambda * (euclidean_distance^2))`. `intercept` and `lambda` are
#'   parameters to be defined -- see [oneimpact::zoi_functions()] for details.
#'   \item If `exp_decay`, the ZoI follows an exponential decay shape: \cr
#'   `intercept * exp(-lambda * euclidean_distance)`. `intercept` and `lambda` are
#'   parameters to be defined -- see [oneimpact::zoi_functions()] for details.
#'   \item If `bartlett`, `linear_decay`, or `tent_decay`, the ZoI follows a
#'   linear decay shape within the ZoI radius (`radius`).
#'   \item If `threshold` or `step`, a constant influence is considered within the
#'   zone of influence radius (`radius`). All pixels closer than
#'   `radius` to infrastructure are considered as "under the influence" of
#'   the nearest feature, with a constant influence value defined by the
#'   `intercept` parameter, and all other pixels are assumed to have
#'   zero influence.
#'   \item If `euclidean`, the function returns the Euclidean distance as a
#'   proxy for the ZoI, even though a proper zone of influence is not defined
#'   in this case.
#'   \item If `log`, the function returns the log-distance:
#'   `log(euclidean_distance, base = log_base)` as a
#'   proxy for the ZoI, even though a proper zone of influence is not defined
#'   in this case.
#'   \item If `sqrt`, the functions returns the square rooted distance:
#'   `sqrt(euclidean_distance)` as a
#'   proxy for the ZoI, even though a proper zone of influence is not defined
#'   in this case.
#' }
#' See details below.
#' Other options still to be implemented (such as other functions and a
#' generic user-defined ZoI function as input).
#'
#' @param where `[character(1)="R"]{"R", "GRASS"}` \cr Where should the
#' computation be done? Default is `"R"`. If `where = "GRASS"`, the R session
#' must be linked to an open GRASS GIS session in a specific location and mapset.
#'
#' @param intercept `[numeric(1)=1]` \cr Maximum value of the ZoI functions at
#' when the distance from disturbance sources is zero (`x = 0`).
#' For the `threshold_decay` and `step_decay` functions, `intercept` is
#' the constant value of the Zone of Influence within the ZoI `radius`.
#' For the other ZoI functions, `intercept`
#' is the value of the functions at the origin (where the sources of disturbance
#' are located, i.e. `x = 0`).
#' Default is `intercept = 1`. This parameter is
#' ignored if `type = "euclidean"`, `type = "log"`, or `type = "sqrt"`.
#'
#' @param zoi_limit `[numeric(1)=0.05]` \cr For non-vanishing functions
#' (e.g. `exp_decay`, `gaussian_decay`), this value is used to set the relationship
#' between the ZoI radius and the decay functions:
#' `radius` is defined as the minimum distance at which the ZoI assumes values
#' below `zoi_limit`. The default is 0.05. This parameter is used only
#' if `radius` is not `NULL`.
#'
#' @param lambda `[numeric(2)=NULL]` \cr For the Gaussian and exponential decay
#' functions (`type = "Gauss"` and `type = "exp_decay"`), `lambda` is the decay
#' parameter of the function. Notice that the interpretation of `lambda` is different depending on the
#' the function -- see [oneimpact::zoi_functions()] for definitions.
#' For the Gaussian decay function, the value for `lambda` is only considered if both
#' `radius = NULL` and `sigma = NULL`. For the exponential decay function,
#' the value for `lambda` is only considered if both `radius = NULL` and `half_life = NULL`.
#
#' @param log_base `[numeric(1)=exp(1)]` \cr Base of the logarithm, if `type = log`.
#'
#' @param dist_offset `[numeric(1)=0]` \cr Number to add to the Euclidean
#' distance before transforming it, to avoid `-Inf`/`Inf` values
#' (e.g. in the case of log transformation). It should be a very small value
#' compared to the range of values of Euclidean distance, not to influence any
#' further analyses.
#'
#' @param zeroAsNA `[logical(1)=FALSE]` \cr If `TRUE` treats cells that are
#' zero as if they were `NA`. Only used for computations in R (`where = R`).
#'
#' @param extent_x_cut,entent_y_cut `[numeric vector(2)=NULL]` \cr Vector
#' representing the minimum and
#' maximum extent in x and y for the final output, in the format c(min,max).
#' It is intended to keep only a region of interest, for standardizing the
#' parameters and region when comparing the resulting ZoI maps with the
#' cumulative ZoI, calculated through [oneimpact::calc_zoi_cumulative()].
#' If `NULL` (default), this parameter is ignored.
#'
#' @param g_output_map_name `[character(1)=NULL]` \cr Name of the output map name,
#' to be used only within GRASS (if `where = "GRASS"`). By default, this is `NULL`
#' and the output map names are a concatenation of the input map name
#' (e.g. `"map_houses"`) and the decay function and radius used
#' (e.g. for `type = "exp_decay"` and `radius = 1000`, the name would be
#' `"map_houses_exp_decay1000"`).
#' This parameter is ignored when the calculations are performed in R
#' (`where = "R"`).
#'
#' @param g_dist_metric `[character(1)="euclidean"]{"euclidean", "geodesic",
#' "squared", "maximum", "manhattan"}` \cr
#' If the calculations are perfomed within GRASS GIS, this is the `metric`
#' argument to calculate the distance
#' from the infrastructure features with the module `r.grow.distance`.
#' More information at the GRASS GIS documentation for
#' [this function](https://grass.osgeo.org/grass82/manuals/r.grow.distance.html).
#' This parameter is ignored when the calculations are performed in R
#' (`where = "R"`).
#'
#' @param g_input_as_region `[logical(1)=FALSE]` \cr Should the input map `x` be
#' used to redefine the working region in GRASS before the ZoI calculation?
#' If `TRUE`, `x` is used to define the region with `g.region`. If `FALSE`,
#' the region previously defined in the GRASS GIS session is used for computation.
#' Default is `FALSE`. This parameter is ignored when the calculations are performed in R
#' (`where = "R"`).
#'
#' @param g_remove_intermediate `[logical(1)=TRUE]` \cr Should the intermediate
#' maps created for computing the output map be excluded in the end of the
#' process? Only used when `where = "GRASS"`.
#'
#' @param g_print_expression `[logical(1)=FALSE]` \cr Should the expression for
#' transforming the raster of distance into ZoI values should be printed in the prompt?
#' Only used when `where = "GRASS"` and `verbose = TRUE` for debugging the
#' result of `r.mapcalc`.
#'
#' @param g_overwrite `[logical(1)=FALSE]` \cr If the a map already exists with the
#' name `g_output_map_name` in the working GRASS GIS location and mapset, should
#' it be overwritten? Only used when `where = "GRASS"`.
#'
#' @param verbose `[logical(1)=FALSE]` \cr Should messages of the computation steps
#' be printed in the prompt along the computation?
#'
#' @param ... \cr Adittional parameters passed to [terra::distance()]
#' or to the ZoI functions (see [oneimpact::zoi_functions()]) when the
#' calculations are performed in R.
#' No additional parameters implemented for computation in GRASS GIS.
#'
#' @returns If the calculations are performed in R (`where = "R"`), the function
#' returns a `RasterLayer` (or `SpatRaster`, according to the class of the input
#' object) with the zone of influence of the nearest feature. If multiple values
#' of `radius` are provided, a stack of rasters is returned.
#' If the calculations are performed in GRASS GIS (`where = "GRASS"`),
#' the maps are kept only within the GRASS GIS location/mapset and the function
#' returns the name of the calculated maps. \cr
#' If the computation is done in GRASS GIS, the output is name of
#' the output raster map(s) within the GRASS GIS location and mapset of the
#' current session. The user can retrieve these maps to R using
#' [rgrass::read_RAST] or export them outside GRASS using the
#' `r.out.gdal` module, for instance.
#'
#' @seealso See [oneimpact::zoi_functions()] for some ZoI function shapes. \cr
#' See also [terra::distance()] for details on the calculation of the distance
#' to the nearest future in R. \cr
#' See
#' [r.grow.distance](https://grass.osgeo.org/grass80/manuals/r.mfilter.html),
#' for details on the calculation of the distance to the nearest future in GRASS. \cr
#' See [oneimpact::calc_zoi_cumulative()] for the computation of the cumulative
#' zone of influence and density of multiple features.
#'
#' @example examples/calc_zoi_nearest_example.R
#' @example examples/calc_zoi_nearest_grass_example.R
#'
#' @export
calc_zoi_nearest <- function(
  x,
  radius = NULL,
  type = c("Gauss", "exp_decay", "bartlett", "half_norm", "threshold", "step",
           "euclidean", "log", "sqrt")[1],
  where = c("R", "GRASS")[1],
  intercept = 1,
  zoi_limit = 0.05,
  lambda = NULL,
  log_base = exp(1),
  dist_offset = 0,
  zeroAsNA = FALSE,
  extent_x_cut = NULL,
  extent_y_cut = NULL,
  g_output_map_name = NULL,
  g_dist_metric = c("euclidean", "geodesic", "squared", "maximum", "manhattan")[1],
  g_input_as_region = FALSE,
  g_remove_intermediate = TRUE,
  g_print_expression = FALSE,
  g_overwrite = FALSE,
  verbose = FALSE,
  ...) {

  # intercept and constant influence could be only one parameter

  # Run in R
  if(where %in% c("R", "r")) {
    if(is.null(extent_x_cut)) extent_x_cut <- terra::ext(x)[c(1,2)]
    if(is.null(extent_y_cut)) extent_y_cut <- terra::ext(x)[c(3,4)]

    zoi_nearest <- calc_zoi_nearest_r(x = x,
                                      radius = radius,
                                      type = type,
                                      intercept = intercept,
                                      zoi_limit = zoi_limit,
                                      lambda = lambda,
                                      log_base = log_base,
                                      dist_offset = dist_offset,
                                      zeroAsNA = zeroAsNA,
                                      extent_x_cut = extent_x_cut,
                                      extent_y_cut = extent_y_cut,
                                      verbose = verbose,
                                      ...)

    return(zoi_nearest)
  } else {

    # Run in GRASS GIS
    if(where %in% c("GRASS", "grass", "GRASS GIS", "grass gis")) {
      zoi_nearest <- calc_zoi_nearest_grass(x = x,
                                            radius = radius,
                                            type = type,
                                            intercept = intercept,
                                            zoi_limit = zoi_limit,
                                            zeroAsNA = zeroAsNA,
                                            lambda = lambda,
                                            log_base = log_base,
                                            dist_offset = dist_offset,
                                            extent_x_cut = extent_x_cut,
                                            extent_y_cut = extent_y_cut,
                                            g_output_map_name = g_output_map_name,
                                            g_dist_metric = g_dist_metric,
                                            g_input_as_region = g_input_as_region,
                                            g_remove_intermediate = g_remove_intermediate,
                                            g_print_expression = g_print_expression,
                                            g_overwrite = g_overwrite,
                                            verbose = verbose,
                                            ...)

      return(zoi_nearest)
    }
  }

}

# implementation in R
calc_zoi_nearest_r <- function(
  x,
  radius = NULL,
  type = c("euclidean", "log", "sqrt", "exp_decay", "bartlett", "Gauss", "half_norm", "threshold", "step")[1],
  intercept = 1,
  zoi_limit = 0.05,
  zoi_hl_ratio = NULL,
  half_life = NULL,
  lambda = NULL,
  sigma = NULL,
  log_base = exp(1),
  dist_offset = 0,
  zeroAsNA = FALSE,
  extent_x_cut = terra::ext(x)[c(1,2)],
  extent_y_cut = terra::ext(x)[c(3,4)],
  verbose = FALSE,
  ...) {

  # check if the input is a terra or raster object
  if(class(x) %in% c("SpatRaster")) {
    use_terra <- TRUE
  } else {
    if(class(x) %in% c("RasterLayer", "RasterBrick", "RasterStack")) {
      use_terra <- FALSE
    } else {
      classes <- c("SpatRaster", "RasterLayer", "RasterBrick", "RasterStack")
      stop(paste0("Please make sure x is an object of one of these classes: ",
                  paste(classes, collapse = ","), "."))
    }
  }

  # check if the transformation is valid
  possible_types <- c("euclidean", "log", "sqrt", "exp_decay", "bartlett",
                      "bartlett_decay", "linear_decay", "tent_decay",
                      "Gauss", "gaussian_decay", "gauss", "gaussian",
                      "half_norm", "threshold", "step")
  if(!(type %in% possible_types))
    stop("You should select an appropriate method ('type' parameter)
         for the calculation the ZoI of the nearest feature.")

  # check if the input raster presents only a single value (1,NA)
  # if so, transform it into a binary map (1,0)
  r0 <- x
  if(zeroAsNA) {
    if(use_terra) {
      r0 <- terra::classify(r0, cbind(0, NA)) # binary map
    } else {
      r0 <- raster::reclassify(r0, cbind(0, NA)) # binary map
    }
  }

  # Euclidean distance
  if(verbose) print(paste0("Computing the Euclidean distance to the nearest feature..."))
  dist_r <- terra::distance(r0, ...)

  # transform distance
  if(type != "euclidean")
    if(type == "log") dist_r <- log(dist_r+dist_offset, base = log_base) else
      if(type == "sqrt") dist_r <- sqrt(dist_r+dist_offset) else
        if(type == "exp_decay") {
          if(verbose) print("Computing the ZoI of the nearest feature with
                            exponential decay shape...")
          dist_r <- exp_decay(x = dist_r,
                              radius = radius,
                              intercept = intercept,
                              lambda = lambda,
                              zoi_limit = zoi_limit,
                              half_life = half_life,
                              zoi_hl_ratio = zoi_hl_ratio)
        } else
          if(type %in% c("bartlett", "Bartlett", "bartlett_decay",
                         "linear_decay", "tent_decay")) {
            if(verbose) print("Computing the ZoI of the nearest feature with
                              linear decay shape...")
            dist_r <- linear_decay(dist_r,
                                   radius = radius,
                                   intercept = intercept)

          } else
            if(type %in% c("Gauss", "half_norm", "gauss",
                           "gaussian_decay", "normal_decay")) {
              if(verbose) print("Computing the ZoI of the nearest feature with
                                Gaussian decay shape...")
              dist_r <- gaussian_decay(x = dist_r,
                                       radius = radius,
                                       intercept = intercept,
                                       lambda = lambda,
                                       zoi_limit = zoi_limit,
                                       sigma = sigma)
            } else {
              if(type %in% c("threshold", "step",
                             "threshold_decay", "step_decay")) {
                if(verbose) print("Computing the ZoI of the nearest feature with
                                  threshold shape...")
                dist_r <- threshold_decay(x = dist_r,
                                          radius = radius,
                                          intercept = intercept)

              } else
                stop("You should select an appropriate transformation method ('type' parameter) for the influence.")
            }

  # rename nearest influence layer, including transformation
  name <- paste0("zoi_nearest", "_", type)
  # and the radius
  zoi_methods <- c("exp_decay", "bartlett", "Gauss",
                   "half_norm", "threshold", "step")
  if(type %in% zoi_methods) name <- paste0(name, radius)
  names(dist_r) <- name

  # # return cropped distance layer
  if(use_terra) {
    terra::crop(dist_r, terra::ext(c(extent_x_cut, extent_y_cut)))
  } else
    raster::crop(dist_r, raster::extent(c(extent_x_cut, extent_y_cut)))
}

# implementation in GRASS
calc_zoi_nearest_grass <- function(
  x,
  radius = NULL,
  type = c("Gauss", "exp_decay", "bartlett", "half_norm", "threshold", "step",
           "euclidean", "log", "sqrt")[1],
  intercept = 1,
  zoi_limit = 0.05,
  zeroAsNA = FALSE,
  zoi_hl_ratio = NULL,
  half_life = NULL,
  lambda = NULL,
  sigma = NULL,
  log_base = exp(1),
  dist_offset = 1,
  extent_x_cut = NULL,
  extent_y_cut = NULL,
  g_output_map_name = NULL,
  g_dist_metric = c("euclidean", "geodesic", "squared", "maximum", "manhattan")[1],
  g_input_as_region = FALSE,
  g_remove_intermediate = TRUE,
  g_print_expression = FALSE,
  verbose = FALSE,
  g_overwrite = FALSE,

  ...) {

  # check if the transformation is valid
  possible_types <- c("euclidean", "log", "sqrt", "exp_decay", "bartlett",
                      "bartlett_decay", "linear_decay", "tent_decay",
                      "Gauss", "gaussian_decay", "gauss", "gaussian",
                      "half_norm", "threshold", "step")
  if(!(type %in% possible_types))
    stop("You should select an appropriate method ('type' parameter)
         for the calculation the ZoI of the nearest feature.")

  # flags
  flags <- c()
  if(!verbose) flags <- c(flags, "quiet")
  if(g_overwrite) flags <- c(flags, "overwrite")

  # flags for g.region
  flags_region <- c("a")
  if(verbose) flags_region <- c(flags_region, "p")

  # intermediate maps to remove
  if(g_remove_intermediate) to_remove <- c()

  # 1. check if there is already a connection with GRASS GIS
  # 2. check if the map is already in GRASS GIS mapset, or if it should be
  # uploaded from the disc or from R
  # check if x is a string that exists within GRASS GIS mapset

  # Start by calculating the Euclidean distance from features

  # output name within GRASS GIS
  out_euclidean <- paste0(x, "_zoi_nearest_euclidean")
  # if there is not transformation
  if(type == "euclidean") {
    # if the user provides an output map name, overwrite it
    if(!is.null(g_output_map_name)) out_euclidean <- g_output_map_name
    # if no transformation, the output influence is the euclidean distance
    out_influence <- out_euclidean
  }

  # region
  if(g_input_as_region)
    rgrass::execGRASS("g.region", raster = x, flags = flags_region)
  if(!is.null(extent_x_cut))
    rgrass::execGRASS("g.region", w = extent_x_cut[1], e = extent_x_cut[2],
                       align = x, flags = flags_region)
  if(!is.null(extent_y_cut))
    rgrass::execGRASS("g.region", s = extent_y_cut[1], n = extent_y_cut[2],
                       align = x, flags = flags_region)

  # check if zero should be replaced by NA
  if(zeroAsNA) {
    # print message
    if(verbose) rgrass::execGRASS("g.message", message = "Setting zeros as NULL data...")
    # copy map and use r.null
    temp_null_map <- "temp_set_null_input"
    # copy map
    rgrass::execGRASS("g.copy", raster = paste0(x, ",", temp_null_map), flags = flags)
    # set nulls
    rgrass::execGRASS("r.null", map = temp_null_map, setnull = "0")
    # check if it should be deleted afterwards
    if(g_remove_intermediate)
      to_remove <- c(to_remove, temp_null_map)
  }

  # print message
  if(verbose) rgrass::execGRASS("g.message", message = "Calculating Euclidean distance...")
  # set input
  if(zeroAsNA) input_euc <- temp_null_map else input_euc <- x
  # calculate
  rgrass::execGRASS("r.grow.distance", input = input_euc, distance = out_euclidean,
                     metric = g_dist_metric, flags = flags)
  # check if it should be deleted afterwards
  # check if the layer already exists first; if not, remove it in the end.
  if(g_remove_intermediate)
    if(type != "euclidean") to_remove <- c(to_remove, out_euclidean)

  # transform distance
  if(type != "euclidean") {

    if(type == "log") {

      # log distance
      # message
      out_message <- "Calculating log-distance..."
      # expression
      expression_influence <- sprintf("log(A + %f, %f)", dist_offset, log_base)

    }

    if(type == "sqrt") {

      # sqrt distance
      # message
      out_message <- "Calculating sqrt-distance..."
      # expression
      expression_influence <- sprintf("sqrt(A + %f)", dist_offset)

    }

    if(type == "exp_decay") {

      # define lambda depending on the input parameter
      # if radius is given, it is used
      if(!is.null(radius)) {
        # if zoi_hl_ratio is null, use zoi_limit
        if(is.null(zoi_hl_ratio)) {
          lambda <- log(1/zoi_limit) / radius
        } else {
          # if zoi_hl_ratio is given, calculate lambda
          half_life <- radius/zoi_hl_ratio
          lambda <- log(2)/half_life
        }

      } else {
        # if radius is not given:
        # and half life is given
        if(!is.null(half_life)) {
          lambda <- log(2)/half_life
        } else {
          # otherwise take it from the parameters
          lambda <- lambda
        }
      }

      # exponential decay influence
      # message
      out_message <- "Calculating exponential decay influence..."
      # expression
      # expression_influence <- sprintf("%f * exp(-%f * A)", intercept, lambda)
      # alternative parameterization with inv_lambda
      inv_lambda <- 1/lambda
      expression_influence <- sprintf("%f * exp(- (1/%f) * A)", intercept, inv_lambda)

    }

    if(type %in% c("bartlett", "Bartlett", "bartlett_decay",
                   "linear_decay", "tent_decay")) {

      # betlett (tent-shaped or linear decay) influence
      # message
      out_message <- "Calculating Bartlett influence..."
      # expression
      expression_influence <- sprintf("if(A <= %f, 1 - (1/%f) * A, 0)", radius, radius)

    }

    if(type %in% c("Gauss", "half_norm", "gauss",
                   "gaussian_decay", "normal_decay")) {
      # define radius or half life, depending on which is given as input

      if(!is.null(radius)) {
        lambda = log(1/zoi_limit) / (radius**2)
      } else {
        if(!is.null(sigma)) {
          lambda = 1/(2*sigma**2)
        } else {
          lambda <- lambda
        }

      }

      # exponential decay influence
      # message
      out_message <- "Calculating half-normal decay influence..."
      # expression
      # expression_influence <- sprintf("%f * exp(-%f * pow(A, 2))", intercept, lambda)
      # alternative parameteization with inv_lambda
      inv_lambda <- 1/lambda
      expression_influence <- sprintf("%f * exp(- (1/%f) * pow(A, 2))", intercept, inv_lambda)

    }

    if(type %in% c("threshold", "step", "threshold_decay", "step_decay")) {

      # threshold influence
      # message
      out_message <- "Calculating threshold influence..."
      # expression
      expression_influence <- sprintf("if(A < %f, %f, 0)", radius, intercept)

    }

    # if the user provides an output map name, use it
    if(!is.null(g_output_map_name))
      out_influence <- g_output_map_name
    else {
      # otherwise, name it according to the method
      out_influence <- paste0(x, "_zoi_nearest_", type)
      # and maybe the radius
      zoi_methods <- c("exp_decay", "bartlett", "Gauss",
                       "half_norm", "threshold", "step")
      if(type %in% zoi_methods) out_influence <- paste0(out_influence, radius)
    }

    # print message
    if(verbose) rgrass::execGRASS("g.message", message = out_message)
    # set region
    if(g_input_as_region)
      rgrass::execGRASS("g.region", raster = out_euclidean, flags = flags_region)

    # compute ZoI
    for(ii in seq_along(expression_influence)) {

      if(g_print_expression) rgrass::execGRASS("g.message", message = expression_influence[ii])
      rgrass::execGRASS("r.mapcalc.simple", expression = expression_influence[ii],
                         a = out_euclidean, output = out_influence[ii], flags = flags)
    }

  }

  # remove intermediate maps
  remove_flags = ifelse(verbose, "f", c("f", "quiet"))
  if(g_remove_intermediate)
    if(length(to_remove) > 0)
      rgrass::execGRASS("g.remove", type = "rast", name = to_remove,
                         flags = remove_flags)

  # return only names
  return(out_influence)
}
