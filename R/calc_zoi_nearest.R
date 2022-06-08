#' Calculate the zone of influence from the nearest feature
#'
#' @description
#' This function takes in a raster with locations or counts of
#' infrastructure and calculates a raster (or set of rasters, in case there is
#' more the one value for `zoi_radius`) representing the zone of influence (ZoI)
#' from the neareast feature of that type of infrastructure. Zones of influence
#' are defined by functions that decay with the Euclidean distance from each
#' infrastructure and their rate of decay is controlled by the ZoI radius
#' (`zoi_radius`), which defines how far the influence of an infrastructure
#' feature goes. By default, the Gaussian decay ZoI is calculated, byt other
#' decay functions might be used (see [oneimpact::zoi_funtions] for examples).
#' The function might also return the Euclidean distance to the nearest feature
#' or a transformation from it (e.g. log- and sqrt-distance from the nearest
#' feature).
#'
#' The procedure might be computed in both R and GRASS GIS. In GRASS, it
#' requires an active connection between the R session and a GRASS GIS
#' location and mapset (through the package [rgrass7]), and that the input
#' maps are already loaded within this GRASS GIS mapset.
#' If the calculations are done in R, the input is a (set of) raster map(s)
#' and the function returns another (set of) raster map(s). If the calculations
#' are done within GRASS GIS, the input is the name of a raster map already
#' loaded in a GRASS GIS location and mapset, and the function returns
#' only the name of the output map. This map is stored in the the GRASS GIS
#' location/mapset, and might be retrieved to R through the
#' [rgrass7::read_RAST] function or exported outside GRASS using the
#' `r.out.gdal` module, for instance.
#'
#' @details
#' In practice, the function first calculated the Euclidean distance from each
#' pixel to the nearest feature and then transforms it according to the ZoI
#' functions. In R, `calc_zoi_nearest` makes use of the [terra::distance]
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
#' preferentially with value 1 (or any other positive value) where the features
#' are located and NA elsewhere. Alternatively, `x` might be a binary (dummy)
#' spatial variable representing the presence of linear or area features, with
#' NA/no-data as background.
#' `x` can be a `RasterLayer` from [raster] package or a [SpatRaster] from
#' [terra] package. If `where = "GRASS"`, `x` must be a string corresponding
#' to the name of the input map within a GRASS GIS location and mapset.
#' Maps without NA as background might be prepared as input for `calc_zoi_nearest`
#' through [raster algebra](https://rspatial.org/terra/pkg/4-algebra.html) in R
#' and e.g. through the use of the module
#' [`r.null`](https://grass.osgeo.org/grass80/manuals/r.null.html) in GRASS GIS.
#'
#' @param zoi_radius `[numeric(1)]` \cr Zone of Influence (ZoI) radius,
#' the distance at which the ZoI vanishes or goes below a given minimum limit value
#' `zoi_limit`. See [oneimpact::zoi_functions] for details. This parameter is
#' ignored if `type = "euclidean"`, `type = "log"`, or `type = "sqrt"`.
#'
#' @param type `[character(1)="Gauss"]{"Gauss", "exp_decay", "bartlett",
#' "threshold", "step", "euclidean", "log","sqrt"}` \cr
#' \itemize{
#'   \item If `Gauss` or `half_norm`, the ZoI follows a half-normal shape: \cr
#'   `N_0 * exp(-lambda * (euclidean_distance^2))`. `N_0` and `lambda` are
#'   parameters to be defined -- see [oneimpact::zoi_functions] for details.
#'   \item If `exp_decay`, the ZoI follows an exponential decay shape: \cr
#'   `N_0 * exp(-lambda * euclidean_distance)`. `N_0` and `lambda` are
#'   parameters to be defined -- see [oneimpact::zoi_functions] for details.
#'   \item If `bartlett`, `linear_decay`, or `tent_decay`, the ZoI follows a
#'   linear decay shape within the ZoI radius (`zoi_radius`).
#'   \item If `threshold` or `step`, a constant influence is consider within the
#'   zone of influence radius (`zoi_radius`). All pixels closer than
#'   `zoi_radius` to infrastructure are considered as "under the influence" of
#'   the nearest feature, with a constant influence value defined by the
#'   `constant_influence` parameter, and all other pixels are assumed to have
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
#' @param log_base `[numeric(1)=exp(1)]` \cr Base of the logarithm, if `type = log`.
#'
#' @param zoi_limit `[numeric(1)=0.05]` \cr For non-vanishing functions
#' (e.g. `exp_decay`, `gaussian_decay`), this value is used to set the relationship
#' between the ZoI radius and the decay functions:
#' `zoi_radius` is defined as the minimum distance at which the ZoI assumes values
#' below `zoi_limit`. The default is 0.05. This parameter is used only
#' if `zoi_radius` is not `NULL`.
#'
#' @param exp_decay_parms `[numeric(2)=c(1,0.01)]` \cr Parameters (`N_0`, `lambda`)
#' for the exponential decay ZoI, if `type = exp_decay`. The value of `lambda`
#' defined here is used only if `zoi_radius = NULL` and `half_life = NULL`,
#' otherwise one of these parameters is used to determine `lambda`.
#' By default, `N_0` is defined as 1, which means the ZoI is 1 where the
#' infrastructure feature is located, and it decreases as the Euclidean distance
#' from it increases.
#'
#' @param hnorm_decay_parms `[numeric(2)=c(1,20)]` \cr Parameters (`N_0`, `sigma`)
#' for the half-normal decay ZoI, if `type = Gauss` or `type = half_normal`.
#' The value of `sigma` defined here is used to define the decay rate `lambda`
#' only if `zoi_radius = NULL` and `half_life = NULL`, otherwise one of these
#' parameters is used to determine `lambda`. By default, `N_0` is defined as 1,
#' which means the ZoI is 1 where the infrastructure feature is located,
#' and it decreases as the Euclidean distance from it increases.
#'
#' @param dist_offset `[numeric(1)=1]` \cr Number to add to the Euclidean
#' distance before transforming it, to avoid `-Inf`/`Inf` values
#' (e.g. in the case of log transformation). It should be a very small value
#' compared to the range of values of Euclidean distance, not to influence any
#' further analyses.
#'
#' @param constant_influence `[numeric(1)=1]` \cr Constant value of the
#' ZoI of the nearest feature if `type = "threshold"` or `type = "step"`.
#' Default is 1. In this case, all pixels closer to any
#' infrastructure than the `zoi_radius` are classified with this constant value.
#'
#' @param extent_x_cut,entent_y_cut `[numeric vector(2)=c(0,1)]` \cr Vector
#' representing the minimum and
#' maximum extent in x and y for the final output, in the format c(min,max).
#' It is intended to keep only a region of interest, for standardizing the
#' parameters and region when comparing the resulting ZoI maps with the
#' cumulative ZoI, calculated through [oneimpact::calc_zoi_cumulative].
#'
#' @param plotit `[logical(1)=FALSE]` \cr Should the outputs be plotted along
#' the calculation? Only used when `where = "R"`.
#'
#' @param output_map_name `[character(1)=NULL]` \cr Name of the output map name,
#' to be used only within GRASS (if `where = "GRASS"`). By default, this is `NULL`
#' and the output map names are a concatenation of the input map name
#' (e.g. "map_houses") and the decay function and zoi_radius used
#' (e.g. for `type = "exp_decay"` and `zoi_radius = 1000`, the name would be
#' "map_houses_exp_decay_1000").
#' This parameter is ignored when the calculations are performed in R
#' (`where = "R"`).
#' @param metric `[character(1)="euclidean"]{"euclidean", "geodesic",
#' "squared", "maximum", "manhattan"}` \cr
#' If the calculations are perfomed within GRASS GIS, this is the `metric`
#' argument to calculate the distance
#' from the infrastructure features with the module `r.grow.distance`.
#' More information at the GRASS GIS documentation for
#' [this function](https://grass.osgeo.org/grass80/manuals/r.grow.distance.html).
#' This parameter is ignored when the calculations are performed in R
#' (`where = "R"`).
#' @param input_as_region `[logical(1)=TRUE]` \cr Should the input map `x` be
#' used to redefine the working GRASS region before cumulative ZoI calculation?
#' If `TRUE`, `x` is used to define the region with `g.region`. If `FALSE`,
#' the region previously defined in the GRASS GIS session is used for computation.
#' @param remove_intermediate `[logical(1)=TRUE]` \cr Should the intermediate
#' maps created for computing the output map be excluded in the end of the
#' process? Only used when `where = "GRASS"`.
#' @param print_expression `[logical(1)=FALSE]` \cr Should the expression for
#' transforming the raster of distance should be printed in the prompt?
#' Only used when `where = "GRASS"` for debugging the result of `r.mapcalc`.
#' @param overwrite `[logical(1)=FALSE]` \cr If the a map already exists with the
#' name `output_map_name` in the working GRASS GIS location and mapset, should
#' it be overwritten? Only used when `where = "GRASS"`.
#' @param quiet `[logical(1)=TRUE]` \cr Should GRASS GIS messages be ommited
#' from the prompt along the computation? Only used when `where = "GRASS"`.
#'
#' @param ... \cr Adittional parameters passed to [terra::distance()]
#' or to the ZoI functions (see [oneimpact::zoi_functions]) when the
#' calculations are performed in R.
#' No additional parameters implemented for computation in GRASS GIS.
#'
#' @returns If the calculations are performed in R (`where = "R"`), the function
#' returns a `RasterLayer` (or `SpatRaster`, according to the class of the input
#' object) with the zone of influence of the nearest feature. If multiple values
#' of `zoi_radius` are providade, a stack of rasters is returned.
#' If the calculations are performed in GRASS GIS (`where = "GRASS"`),
#' the maps are kept only within the GRASS GIS location/mapset and the function
#' returns the name of the calculated maps. \cr
#' If the computation is done in GRASS GIS, the output is name of
#' the output raster map within the GRASS GIS location and mapset of the
#' current session. The user can retrieve these maps to R using
#' [rgrass7::read_RAST] or export them outside GRASS using the
#' `r.out.gdal` module, for instance.
#'
#' @seealso See [oneimpact::zoi_functions] for some ZoI function shapes. \cr
#' See also [terra::distance] for details on the calculation of the distance
#' to the nearest future in R. \cr
#' See
#' [r.grow.distance](https://grass.osgeo.org/grass80/manuals/r.mfilter.html),
#' for details on the calculation of the distance to the nearest future in GRASS. \cr
#' See [oneimpact::calc_zoi_cumulative] for the computation of the cumulative
#' zone of influence and density of multiple features.
#'
#' @example examples/calc_zoi_nearest_example.R
#' @example examples/calc_zoi_nearest_grass_example.R
#'
#' @export
calc_zoi_nearest <- function(
  x,
  zoi_radius = NULL,
  type = c("Gauss", "exp_decay", "bartlett", "half_norm", "threshold", "step",
           "euclidean", "log", "sqrt")[1],
  where = c("R", "GRASS")[1],
  log_base = exp(1),
  zoi_limit = 0.05,
  zoi_hl_ratio = NULL,
  half_life = NULL,
  exp_decay_parms = c(1, 0.01),
  hnorm_decay_parms = c(1, 20),
  constant_influence = 1,
  dist_offset = 0,
  extent_x_cut = NULL,
  extent_y_cut = NULL,
  plotit = FALSE,
  output_map_name = NULL,
  metric = c("euclidean", "geodesic", "squared", "maximum", "manhattan")[1],
  input_as_region = TRUE,
  remove_intermediate = TRUE,
  print_expression = TRUE,
  quiet = TRUE, overwrite = FALSE,
  ...) {

  # Run in R
  if(where %in% c("R", "r")) {
    if(is.null(extent_x_cut)) extent_x_cut <- terra::ext(x)[c(1,2)]
    if(is.null(extent_y_cut)) extent_y_cut <- terra::ext(x)[c(3,4)]

    zoi_nearest <- calc_zoi_nearest_r(x = x,
                                      zoi_radius = zoi_radius,
                                      type = type,
                                      log_base = log_base,
                                      zoi_limit = zoi_limit,
                                      exp_decay_parms = exp_decay_parms,
                                      hnorm_decay_parms = hnorm_decay_parms,
                                      sigma = sigma,
                                      constant_influence = constant_influence,
                                      dist_offset = dist_offset,
                                      extent_x_cut = extent_x_cut,
                                      extent_y_cut = extent_y_cut,
                                      plotit = plotit, ...)

    return(zoi_nearest)
  } else {

    # Run in GRASS GIS
    if(where %in% c("GRASS", "grass", "GRASS GIS", "grass gis")) {
      zoi_nearest <- calc_zoi_nearest_grass(x = x,
                                            zoi_radius = zoi_radius,
                                            type = type,
                                            log_base = log_base,
                                            zoi_limit = zoi_limit,
                                            exp_decay_parms = exp_decay_parms,
                                            hnorm_decay_parms = hnorm_decay_parms,
                                            constant_influence = constant_influence,
                                            dist_offset = dist_offset,
                                            extent_x_cut = extent_x_cut,
                                            extent_y_cut = extent_y_cut,
                                            output_map_name = output_map_name,
                                            metric = metric,
                                            input_as_region = input_as_region,
                                            remove_intermediate = remove_intermediate,
                                            print_expression = print_expression,
                                            quiet = quiet,
                                            overwrite = overwrite,
                                            ...)

      return(zoi_nearest)
    }
  }

}

# implementation in R
calc_zoi_nearest_r <- function(
  x,
  zoi_radius = NULL,
  type = c("euclidean", "log", "sqrt", "exp_decay", "bartlett", "Gauss", "half_norm", "threshold", "step")[1],
  log_base = exp(1),
  zoi_limit = 0.05,
  zoi_hl_ratio = NULL,
  half_life = NULL,
  exp_decay_parms = c(1, 0.01),
  hnorm_decay_parms = c(1, 20),
  sigma = NULL,
  constant_influence = 1,
  dist_offset = 0,
  extent_x_cut = terra::ext(x)[c(1,2)],
  extent_y_cut = terra::ext(x)[c(3,4)],
  plotit = FALSE, ...) {

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

  # Euclidean distance
  dist_r <- terra::distance(x, ...)

  # transform distance
  if(type != "euclidean")
    if(type == "log") dist_r <- log(dist_r+dist_offset, base = log_base) else
      if(type == "sqrt") dist_r <- sqrt(dist_r+dist_offset) else
        if(type == "exp_decay") {
          dist_r <- exp_decay(x = dist_r,
                              zoi_radius = zoi_radius,
                              exp_decay_parms = exp_decay_parms,
                              zoi_limit = zoi_limit,
                              half_life = half_life,
                              zoi_hl_ratio = zoi_hl_ratio)
        } else
          if(type %in% c("bartlett", "Bartlett", "bartlett_decay",
                         "linear_decay", "tent_decay")) {

            ################# change to function later on
            dist_r <- (1 - (1/zoi_radius)*dist_r)
            zero <- dist_r
            size_landscape <- prod(dim(zero)[c(1:2)])
            values(zero) <- rep(0, size_landscape)
            dist_zero_stack <- c(dist_r, zero)

            if(use_terra) {
              dist_r <- terra::app(dist_zero_stack, "max")
            } else {
              dist_zero_stack <- raster::stack(dist_zero_stack)
              dist_r <- raster::calc(dist_zero_stack, max)
            }

          } else
            if(type %in% c("Gauss", "half_norm", "gauss",
                           "gaussian_decay", "normal_decay")) {
              dist_r <- gaussian_decay(x = dist_r,
                                       zoi_radius = zoi_radius,
                                       hnorm_decay_parms = hnorm_decay_parms,
                                       zoi_limit = zoi_limit,
                                       sigma = sigma)
            } else {
              if(type %in% c("threshold", "step",
                             "threshold_decay", "step_decay")) {

                ################# change to function later on
                # for threshold influence of the nearest, keep pixels with
                # dist <= zoi_radius constant
                values(dist_r) <- ifelse(values(dist_r) < zoi_radius,
                                         constant_influence, 0)

              } else
                stop("You should select an appropriate transformation method ('type' parameter) for the influence.")
            }

  # rename nearest influence layer, including transformation
  name <- paste0("zoi_nearest", "_", type)
  # and the zoi_radius
  zoi_methods <- c("exp_decay", "bartlett", "Gauss",
                   "half_norm", "threshold", "step")
  if(type %in% zoi_methods) name <- paste0(name, zoi_radius)
  names(dist_r) <- name

  # should the result be plotted?
  if(plotit) plot(dist_r)

  # # return cropped distance layer
  if(use_terra) {
    terra::crop(dist_r, terra::ext(c(extent_x_cut, extent_y_cut)))
  } else
    raster::crop(dist_r, raster::extent(c(extent_x_cut, extent_y_cut)))
}

# implementation in GRASS
calc_zoi_nearest_grass <- function(
  x,
  zoi_radius = NULL,
  type = c("Gauss", "exp_decay", "bartlett", "half_norm", "threshold", "step",
           "euclidean", "log", "sqrt")[1],
  log_base = exp(1),
  zoi_limit = 0.05,
  zoi_hl_ratio = NULL,
  half_life = NULL,
  exp_decay_parms = c(1, 0.01),
  hnorm_decay_parms = c(1, 20),
  sigma = NULL,
  constant_influence = 1,
  dist_offset = 1,
  extent_x_cut = NULL,
  extent_y_cut = NULL,
  output_map_name = NULL,
  metric = c("euclidean", "geodesic", "squared", "maximum", "manhattan")[1],
  input_as_region = TRUE,
  remove_intermediate = TRUE,
  print_expression = FALSE,
  quiet = FALSE, overwrite = FALSE,
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
  if(quiet) flags <- c(flags, "quiet")
  if(overwrite) flags <- c(flags, "overwrite")

  # flags for g.region
  flags_region <- c("a")
  if(!quiet) flags_region <- c(flags_region, "p")

  # intermediate maps to remove
  if(remove_intermediate) to_remove <- c()

  # 1. check if there is already a connection with GRASS GIS
  # 2. check if the map is already in GRASS GIS mapset, or if it should be
  # uploaded from the disc or from R
  # check if x is a string that exists within GRASS GIS mapset

  # Start by calculating the Euclidean distance from features

  # output name within GRASS GIS
  out_euclidean <- paste0(x, "_inf_nearest_euclidean")
  # if there is not transformation
  if(type == "euclidean") {
    # if the user provides an output map name, overwrite it
    if(!is.null(output_map_name)) out_euclidean <- output_map_name
    # if no transformation, the output influence is the euclidean distance
    out_influence <- out_euclidean
  }

  # region
  if(input_as_region)
    rgrass7::execGRASS("g.region", raster = x, flags = flags_region)
  if(!is.null(extent_x_cut))
    rgrass7::execGRASS("g.region", w = extent_x_cut[1], e = extent_x_cut[2],
                       align = x, flags = flags_region)
  if(!is.null(extent_y_cut))
    rgrass7::execGRASS("g.region", s = extent_y_cut[1], n = extent_y_cut[2],
                       align = x, flags = flags_region)

  # print message
  if(!quiet) rgrass7::execGRASS("g.message", message = "Calculating Euclidean distance...")
  # calculate
  rgrass7::execGRASS("r.grow.distance", input = x, distance = out_euclidean,
                     metric = metric, flags = flags)
  # check if it should be deleted afterwards
  # check if the layer already exists first; if not, remove it in the end.
  if(remove_intermediate)
    if(type != "euclidean") to_remove <- c(to_remove, out_euclidean)

  # transform distance
  if(type != "euclidean") {

    if(type == "log") {

      # log distance
      # message
      out_message <- "Calculating log-distance..."
      # expression
      expression_influence <- sprintf("log(A + %f, %f)", dist_offset, log_base)
      if(print_expression) print(expression_influence)
    }

    if(type == "sqrt") {

      # sqrt distance
      # message
      out_message <- "Calculating sqrt-distance..."
      # expression
      expression_influence <- sprintf("sqrt(A + %f)", dist_offset)
      if(print_expression) print(expression_influence)

    }

    if(type == "exp_decay") {

      # define lambda depending on the input parameter
      # if zoi_radius is given, it is used
      if(!is.null(zoi_radius)) {
        # if zoi_hl_ratio is null, use zoi_limit
        if(is.null(zoi_hl_ratio)) {
          lambda <- log(1/zoi_limit) / zoi_radius
        } else {
          # if zoi_hl_ratio is given, calculate lambda
          half_life <- zoi_radius/zoi_hl_ratio
          lambda <- log(2)/half_life
        }

      } else {
        # if zoi_radius is not given:
        # and half life is given
        if(!is.null(half_life)) {
          lambda <- log(2)/half_life
        } else {
          # otherwise take it from the parameters
          lambda <- exp_decay_parms[2]
        }
      }

      # exponential decay influence
      # message
      out_message <- "Calculating exponential decay influence..."
      # expression
      # expression_influence <- sprintf("%f * exp(-%f * A)", exp_decay_parms[1], lambda)
      # alternative parameterization with inv_lambda
      inv_lambda <- 1/lambda
      expression_influence <- sprintf("%f * exp(- (1/%f) * A)", exp_decay_parms[1], inv_lambda)
      if(print_expression) print(expression_influence)

    }

    if(type %in% c("bartlett", "Bartlett", "bartlett_decay",
                    "linear_decay", "tent_decay")) {

      # betlett (tent-shaped or linear decay) influence
      # message
      out_message <- "Calculating Bartlett influence..."
      # expression
      expression_influence <- sprintf("if(A <= %f, 1 - (1/%f) * A, 0)", zoi_radius, zoi_radius)
      if(print_expression) print(expression_influence)

    }

    if(type %in% c("Gauss", "half_norm", "gauss",
                   "gaussian_decay", "normal_decay")) {
      # define zoi_radius or half life, depending on which is given as input

      if(!is.null(zoi_radius)) {
        lambda = log(1/zoi_limit) / (zoi_radius**2)
      } else {
        if(!is.null(sigma)) {
          lambda = 1/(2*sigma**2)
        } else {
          lambda <- hnorm_decay_parms[2]
        }

      }

      # exponential decay influence
      # message
      out_message <- "Calculating half-normal decay influence..."
      # expression
      # expression_influence <- sprintf("%f * exp(-%f * pow(A, 2))", hnorm_decay_parms[1], lambda)
      # alternative parameteization with inv_lambda
      inv_lambda <- 1/lambda
      expression_influence <- sprintf("%f * exp(- (1/%f) * pow(A, 2))", hnorm_decay_parms[1], inv_lambda)
      if(print_expression) print(expression_influence)

    }

    if(type %in% c("threshold", "step", "threshold_decay", "step_decay")) {

      # threshold influence
      # message
      out_message <- "Calculating threshold influence..."
      # expression
      expression_influence <- sprintf("if(A < %f, %f, 0)", zoi_radius, constant_influence)
      if(print_expression) print(expression_influence)
    }

    # if the user provides an output map name, use it
    if(!is.null(output_map_name))
      out_influence <- output_map_name
    else {
      # otherwise, name it according to the method
      out_influence <- paste0(x, "_zoi_nearest_", type)
      # and maybe the zoi_radius
      zoi_methods <- c("exp_decay", "bartlett", "Gauss",
                       "half_norm", "threshold", "step")
      if(type %in% zoi_methods) out_influence <- paste0(out_influence, zoi_radius)
    }

    # print message
    if(!quiet) rgrass7::execGRASS("g.message", message = out_message)
    # set region
    if(input_as_region)
      rgrass7::execGRASS("g.region", raster = out_euclidean, flags = flags_region)
    # calculate
    rgrass7::execGRASS("r.mapcalc.simple", expression = expression_influence,
                       a = out_euclidean, output = out_influence, flags = flags)
  }

  # remove intermediate maps
  remove_flags = ifelse(quiet, c("f", "quiet"), "f")
  if(remove_intermediate & length(to_remove) > 1)
    rgrass7::execGRASS("g.remove", type = "rast", name = to_remove,
                       flags = remove_flags)

  # return only names
  return(out_influence)
}
