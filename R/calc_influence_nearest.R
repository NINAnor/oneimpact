#' Calculate the influence from the nearest feature
#'
#' This function takes in a raster with locations of infrastructure and calculates
#' a raster representing the influence to each pixel from the neareast feature of that type
#' of infrastructure. By default, the output measure of influence is the Euclidean distance
#' to the nearest feature. However, this influence measure can be changed to
#' transformed distance measures (so far, log- and sqrt-distance) or to decay distance functions
#' (so far, exponential decay distance, half-normal distance, and Bartlett distance). For the decay distance
#' measures of influence, a zone of influence (zoi) might be used as input for designing
#' the curve. The calculations might be done within R or within GRASS GIS. In the latter case,
#' it requires an active connection between the R session and a GRASS GIS location and mapset
#' (through the package [rgrass7]), and that the input maps are already loaded within this GRASS GIS mapset.
#'
#' The input raster should have positive values in the pixels where infrastructure
#' are located and NA/no-data in all other places. The input raster is supposed to
#' represent the location of point, line, or polygon infrastructure (e.g. houses, roads, mining areas),
#' but any landscape variable whose representation might be one of those would fit here
#' (e.g. areas of forest or any other habitat type or land cover). We recommend that the
#' input raster has a metric projection, so that distances and influence measures are
#' based on distance to infrastructure measured in meters.
#'
#' If the calculations are done in R, the input is a raster map and the function returns
#' another raster map. If the calculations are done within GRASS GIS, the input is the name
#' of a raster map already loaded in a GRASS GIS location and mapset, and the function returns
#' only the name of the output map. This map is stored in the the GRASS GIS location/mapset,
#' and might be retrieved to R through the [rgrass7::readRAST()] function, for instance.
#' The calculations in GRASS GIS are done through the modules
#' [`r.grow.distance`](https://grass.osgeo.org/grass78/manuals/r.grow.distance.html) for computing
#' the Euclidean distance to the nearest features and
#' [`r.mapcalc.simple`](https://grass.osgeo.org/grass78/manuals/r.mapcalc.simple.html)
#' to transform this distance into the different forms of measures of influence of the nearest features.
#'
#' TO IMPROVE: implement any generic function as input to "transform".
#'
#' @param x `[RasterLayer,SpatRaster]` \cr Raster representing locations of features, with value 1
#' (or any other positive value) where the features are located and NA elsewhere.
#' Can be a [Raster-class] object (e.g. `RasterLayer`) from [raster-package] or a
#' [SpatRaster] from [terra-package].
#' The features might be represented by points (e.g. houses, cabins, wind turbines), lines
#' (e.g. roads, power lines, railways), or polygons (e.g. mining areas, urban areas, or
#' any other spatial variable represented as polygons or areas).
#'
#' @param zoi `[numeric(1)=NULL]` \cr Zone of Influence (ZoI), in map units (preferentially meters).
#' The ZoI is the distance, scale, or buffer size around a feature up to which we consider there is
#' an effect or influence of an infrastructure or variable. It is considered only when
#' `transform = "bartlett"` or `transform = "exp_decay"`. \cr
#' For the Bartlett influence, it corresponds to the distance beyond which the distance is zero. \cr
#' For the exponential decay influence, the ZoI is used to define the half-life and the lambda
#' of the exponential decay function, based on the parameter `zoi_hl_ratio`,
#' that defines the ratio between the ZoI and the half-life. Since the half-life is the value
#' where the exponential decay decreases by `0.5`, a ratio of, for instance, `zoi_hl_ratio = 4` (default)
#' would mean that the ZoI is defined as the value where the exponential decay decreases to `0.5^4 = 0.0625`.
#' In this case, if `zoi = 4000` m, this means that the ZoI is four times higher than the half-life, i.e.
#' `half_life = 1000` and `lambda = log(2)/half_life = 6.93e-4`. The definition of a zone of
#' influence does not imply a cuttoff of the exponential decay function but is only used to define
#' its parameters, based on the defined `zoi_hl_ratio` parameter.
#'
#' @param transform `[character(1)=NULL]{"log","sqrt", "exp_decay", "Gauss", "bartlett"}` \cr
#' By default, NULL, when the measure of influence is the Euclidean distance to the nearest
#' feature.
#' \itemize{
#'   \item If `log`, the influence measure is the log-distance: `log(euclidean_distance, base = log_base)`.
#'   \item If `sqrt`, the influence measure is the square rooted distance: `sqrt(euclidean_distance)`.
#'   \item If `exp_decay`, the influence measure is the exponential decay distance: \cr
#' `N_0 * exp(-lambda * euclidean_distance)`. `N_0` and `lambda` are parameters to be defined.
#' The decayment rate lambda might be defined in terms of the exponential half-life or
#' a definition of zone of influence (ZoI).
#'   \item If `Gauss` or `half_norm`, the influence measure is the half-normal decay distance: \cr
#' `N_0 * exp(-lambda * (euclidean_distance^2))`. `N_0` and `lambda` are parameters to be defined.
#' The decayment rate lambda might be defined in terms of the Gaussian half-life or
#' a definition of zone of influence (ZoI), or yet by the standard deviation (sigma) of the
#' half-normal distribution (`lambda = 0.5/(sigma^2)`).
#'   \item If `bartlett`, the influence measure is a triangular tent-shaped decay distance is returned.
#'   \item If `threshold` or `step`, a constant influence is consider within the zone of influence (ZoI).
#' All pixels closer than `zoi` to infrastructure are considered as "under the influence" of the nearest
#' feature, with a constant influence value defined by the `constant_influence` parameter, and all other
#' pixels are assumed to have zero influence.
#' }
#' See details below.
#' Other options still to be implemented (such as other functions and a generic user-defined
#' function as input).
#'
#' @param where `[character(1)="R"]{"R", "GRASS"}` \cr Where the calculations should be performed, in R or
#' GRASS GIS. Calculations within GRASS GIS (`where = "GRASS"`) require an active connection between the R session
#' and a GRASS GIS location and mapset and that the input map is already loaded in this the current mapset.
#'
#' @param log_base `[numeric(1)=exp(1)]` \cr Base of the logarithm, if `transform = log`.
#'
#' @param zoi_hl_ratio `[numeric(1)=4]` \cr REVIEW THAT, THE IDEA IS SIMILAR BUT SLIGHTLY DIFFERENT FOR GAUSSIAN
#' Ratio between the zone of influence (ZoI) of the infrastrucure
#' and the half-life of the exponential decay distance function. It is used to define the
#' the half-life and decay parameter lambda of the exponential decay function, based on a value of ZoI.
#' By changing this parameter, one can be more or less flexible in the definition of ZoI for
#' the exponential decay influence.
#' For example: if `zoi_hl_ratio = 4`, the ZoI is defined as the distance from the nearest infrastructure
#' where the decay distance influence decreases to `0.5^4 = 0.0625`. If `zoi_hl_ratio = 6`, the ZoI
#' is defined as tyhe distance when the exponential decay influence decreases to `0.5^6 = 0.015625`
#'
#' @param half_life `[numeric(1)=NULL]` \cr Half-life of the exponential decay function or half-normal decay function,
#' in case `transform = exp_decay` or `transform = half_norm`. The exponential decay exponent (lambda) from the exponential function is defined as
#' `lambda = log(2)/half_life`. For the half-normal decay the exponent (lambda) is `lambda = log(2)/(half_life**2)`
#' By definition, when one gets away from the source feature by a
#' distance interval equals to `half_life`, the magnitude of the exponential decay distance decreases by 1/2.
#' This means that, for instance, at a distance of `4*half_life` to the nearest feature, the exponential decay
#' influence has a magnitude of `(1/2)^4 = 1/16 ~ 0.06`. This can be useful to define the
#' Zone of Influence (ZoI) for exponential decay distances. If the `zoi` parameter is not NULL, the
#' `half_life` parameter is ignored and half-life of the exponential decay is defined by the `zoi` and
#' `zoi_hl_ratio` parameters.
#' The half-life definition for the half-normal/Gaussian function is slightly different. While the definition¨
#' of half-life still holds, the relationship with its decay is different. The distance at which
#' the function decreases to `(1/2)^4 = 1/16 ~ 0.06`, for instance, is `sqrt(4)*half_life`.
#'
#' @param exp_decay_parms `[numeric(2)=c(1,0.01)]` \cr Parameters (`N_0`, `lambda`) for the exponential decay
#' influence, if `transform = exp_decay`. The value of `lambda` defined here is used only if
#' `zoi = NULL` and `half_life = NULL`, otherwise one of these parameters is used to determine `lambda`.
#' By default, `N_0` is defined as 1, which means the influence is 1 where the infrastructure feature
#' is located, and it decreases as the Euclidean distance from it increases.
#'
#' @param hnorm_decay_parms `[numeric(2)=c(1,20)]` \cr Parameters (`N_0`, `sigma`) for the half-normal decay
#' influence, if `transform = Gauss` or `transform = half_normal`. The value of `sigme` defined here is used
#' to define the decay rate `lambda`only if `zoi = NULL` and `half_life = NULL`, otherwise one of these
#' parameters is used to determine `lambda`. By default, `N_0` is defined as 1, which means the influence is 1
#' where the infrastructure feature is located, and it decreases as the Euclidean distance from it increases.
#'
#' @param dist_offset `[numeric(1)=1]` \cr Number to add to the Euclidean distance before transforming it,
#' to avoid `-Inf`/`Inf` values (e.g. in the case of log). It should be a very small value compared to the
#' range of values of Euclidean distance, not to influence any further analyses.
#'
#' @param constant_influence `[numeric(1)=1]` \cr Constant value of the influence of the nearest feature if
#' `transform = "threshold"` or `transform = "step"`. Default is 1. In this case, all pixels closer to any
#' infrastructure than the `zoi` are classified with this constant value.
#'
#' @param extent_x_cut,entent_y_cut `[numeric vector(2)=c(0,1)]` \cr Vectors representing the minimum and
#' maximum extent in x and y for the final output, in the format c(min,max). Used to cut the raster to
#' specific smaller extents. The default is to keep the same extent of the input raster.
#' Not used for computation in GRASS GIS yet, but it might be implemented if necessary.
#'
#' @param plotit `[logical(1)=FALSE]` \cr Should the outputs be plotted along the calculation? Valid only when
#' the calculations are done within R (not in GRASS GIS).
#'
#' @param output_map_name `[character(1)=NULL]` \cr Name of the output map name, to be used only within
#' GRASS (if `where = "GRASS"`). Bu default, this is NULL and the output map names are a concatenation of
#' the input map name (e.g. "map_houses") and the decay function and zoi used (e.g. for `transform = "exp_decay"`
#' and `zoi = 1000`, the name would be "map_houses_exp_decay_1000").
#' This parameter is ignored when the calculations are performed in R (if `where = "R"`).
#'
#' @param metric `[character(1)="euclidean"]{"euclidean", "geodesic", "squared", "maximum", "manhattan"}` \cr
#' If the calculations are perfomed within GRASS GIS, this is the `metric` argument to calculate the distance
#' from the infrastructure features with the module `r.grow.distance`. More information at the
#' [GRASS GIS documentation for this function](https://grass.osgeo.org/grass78/manuals/r.grow.distance.html).
#' This parameter is ignored when the calculations are performed in R (if `where = "R"`).
#'
#' @param remove_intermediate = TRUE
#' This parameter is ignored when the calculations are performed in R (if `where = "R"`).
#'
#' @param print_expression = TRUE
#' This parameter is ignored when the calculations are performed in R (if `where = "R"`).
#'
#' @param quiet,overwrite `[logical(1)]` \cr Whether the procedures is GRASS GIS should be run
#' quetly (flag `quiet = TRUE`) and whether the output maps should be overwriten (flag `overwrite = TRUE`).
#' These parameters is ignored when the calculations are performed in R (if `where = "R"`).
#'
#' @param ... \cr Adittional parameters passed to [terra::distance()] when the calculations are performed in R.
#' No additional parameters implemented for computation in GRASS GIS.
#'
#' @returns If the calculations are performed in R (if `where = "R"`), the function returns a `RasterLayer`
#' (or `SpatRaster`, according to the class of the input object) with the
#' influence of the nearest feature. If the calculations are performed in GRASS GIS (if `where = "GRASS"`),
#' the maps are kept only within the GRASS GIS location/mapset and the function returns the name of the
#' calculated maps. \cr
#' By default, the output influence map is the Euclidean distance to the nearest feature.
#' Depending on the choice of `transform`, the output distance can be log- or sqrt-transformed distance,
#' or one can choose to calculate the exponential decay, half-normal decay, or Bartlett decay influence.
#' Other types of transformation to be implemented in the future.
#'
#' @example examples/calc_inf_nearest_example.R
#' @example examples/calc_inf_nearestGRASS_example.R
#'
#' @export
calc_influence_nearest <- function(
  x,
  zoi = NULL,
  transform = NULL, #c("log", "sqrt", "exp_decay", "bartlett", "Gauss", "half_norm", "threshold", "step")[1],
  where = c("R", "GRASS")[1],
  log_base = exp(1),
  zoi_hl_ratio = 4,
  half_life = NULL,
  exp_decay_parms = c(1, 0.01),
  hnorm_decay_parms = c(1, 20),
  constant_influence = 1,
  dist_offset = 1,
  extent_x_cut = NULL,
  extent_y_cut = NULL,
  plotit = FALSE,
  output_map_name = NULL,
  metric = c("euclidean", "geodesic", "squared", "maximum", "manhattan")[1],
  remove_intermediate = TRUE,
  print_expression = TRUE,
  quiet = TRUE, overwrite = FALSE,
  ...) {

  # Run in R
  if(where %in% c("R", "r")) {
    if(is.null(extent_x_cut)) extent_x_cut <- terra::ext(x)[c(1,2)]
    if(is.null(extent_y_cut)) extent_y_cut <- terra::ext(x)[c(3,4)]

    inf_nearest <- calc_influence_nearest_r(x = x,
                                            zoi = zoi,
                                            transform = transform,
                                            log_base = log_base,
                                            zoi_hl_ratio = zoi_hl_ratio,
                                            half_life = half_life,
                                            exp_decay_parms = exp_decay_parms,
                                            hnorm_decay_parms = hnorm_decay_parms,
                                            constant_influence = constant_influence,
                                            dist_offset = dist_offset,
                                            extent_x_cut = extent_x_cut,
                                            extent_y_cut = extent_y_cut,
                                            plotit = plotit, ...)

    return(inf_nearest)
  } else {

    # Run in GRASS GIS
    if(where %in% c("GRASS", "grass", "GRASS GIS", "grass gis")) {
      inf_nearest <- calc_influence_nearest_GRASS(x = x,
                                                  zoi = zoi,
                                                  transform = transform,
                                                  log_base = log_base,
                                                  zoi_hl_ratio = zoi_hl_ratio,
                                                  half_life = half_life,
                                                  exp_decay_parms = exp_decay_parms,
                                                  hnorm_decay_parms = hnorm_decay_parms,
                                                  constant_influence = constant_influence,
                                                  dist_offset = dist_offset,
                                                  extent_x_cut = extent_x_cut,
                                                  extent_y_cut = extent_y_cut,
                                                  output_map_name = output_map_name,
                                                  metric = metric,
                                                  remove_intermediate = remove_intermediate,
                                                  print_expression = print_expression,
                                                  quiet = quiet,
                                                  overwrite = overwrite,
                                                  ...)

      return(inf_nearest)
    }
  }

}

calc_influence_nearest_r <- function(
  x,
  zoi = NULL,
  transform = NULL, #c("log", "sqrt", "exp_decay", "bartlett", "Gauss", "half_norm", "threshold", "step")[1],
  log_base = exp(1),
  zoi_hl_ratio = 4,
  half_life = NULL,
  exp_decay_parms = c(1, 0.01),
  hnorm_decay_parms = c(1, 20),
  constant_influence = 1,
  dist_offset = 1,
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

  # Euclidean distance
  dist_r <- terra::distance(x, ...)

  # transform distance
  if(!is.null(transform))
    if(transform == "log") dist_r <- log(dist_r+dist_offset, base = log_base) else
      if(transform == "sqrt") dist_r <- sqrt(dist_r+dist_offset) else
        if(transform == "exp_decay") {
          # define zoi or half life, depending on which is given as input

          # if zoi is given, it is used.
          # if zoi is not given:
          if(is.null(zoi)) {
            # and half_life is not given either:
            if(is.null(half_life)) {
              # lambda is calculated from exp_decay_parms
              lambda <- exp_decay_parms[2]
            } else {
              # if half_life is given, lambda is calculated from that
              zoi <- half_life * zoi_hl_ratio # not used!
              lambda <- log(2)/half_life
            }
          } else {
            # if zoi is given, it is used to
            # define half_life and lambda
            half_life <- zoi/zoi_hl_ratio
            lambda <- log(2)/half_life
          }

          dist_r <- exp_decay_parms[1] * exp(-lambda * dist_r)
        } else
          if(transform == "bartlett") {
            dist_r <- (1 - (1/zoi)*dist_r)
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
            if(transform %in% c("Gauss", "half_norm")) {
              # define zoi or half life, depending on which is given as input

              # if zoi is given, it is used.
              # if zoi is not given:
              if(is.null(zoi)) {
                # and half_life is not given either:
                if(is.null(half_life)) {
                  # lambda is calculated from exp_decay_parms
                  lambda <- 0.5/(hnorm_decay_parms[2]**2) # lambda = 0.5/(sigma^2)
                } else {
                  # if half_life is given, lambda is calculated from that
                  zoi <- half_life * sqrt(zoi_hl_ratio) # not used!
                  lambda <- log(2)/(half_life**2)
                }
              } else {
                # if zoi is given, it is used to
                # define half_life and lambda
                half_life <- zoi/sqrt(zoi_hl_ratio)
                lambda <- log(2)/(half_life**2)
              }

              dist_r <- hnorm_decay_parms[1] * exp(-lambda * (dist_r**2))
            } else {
              if(transform %in% c("threshold", "step")) {

                # for threshold influence of the nearest, keep pixels with dist <= zoi constant
                values(dist_r) <- ifelse(values(dist_r) < zoi, constant_influence, 0)

              } else
                stop("You should select an appropriate transformation method for the influence.")
            }

  # rename nearest influence layer, including transformation
  name <- ifelse(is.null(transform), "influence_nearest",
                 paste0("influence_nearest", "_", transform))
  # and the zoi
  zoi_methods <- c("exp_decay", "bartlett", "Gauss", "half_norm")
  if(!is.null(transform))
    if(transform %in% zoi_methods) name <- paste0(name, zoi)
  names(dist_r) <- name
  # should the result be plotted?
  if(plotit) plot(dist_r)

  # # return cropped distance layer
  if(use_terra) {
    terra::crop(dist_r, terra::ext(c(extent_x_cut, extent_y_cut)))
  } else
    raster::crop(dist_r, raster::extent(c(extent_x_cut, extent_y_cut)))
}

### clean maps
calc_influence_nearest_GRASS <- function(
  x,
  zoi = NULL,
  transform = NULL,
  log_base = exp(1),
  zoi_hl_ratio = 4,
  half_life = NULL,
  exp_decay_parms = c(1, 0.01),
  hnorm_decay_parms = c(1, 20),
  constant_influence = 1,
  dist_offset = 1,
  extent_x_cut = NULL,
  extent_y_cut = NULL,
  output_map_name = NULL,
  metric = c("euclidean", "geodesic", "squared", "maximum", "manhattan")[1],
  remove_intermediate = TRUE,
  print_expression = FALSE,
  quiet = FALSE, overwrite = FALSE,
  ...) {

  # check if the transformation is valid
  possible_transformations <- c("log", "sqrt", "exp_decay", "bartlett", "Gauss", "half_norm", "threshold", "step")
  if(!is.null(transform))
    if(!(transform %in% possible_transformations))
      stop("You should select an appropriate transformation method for distance.")

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
  # 2. check if the map is already in GRASS GIS mapset, or if it should be uploaded from the disc or from R
  # check if x is a string that exists within GRASS GIS mapset

  # Start by calculating the Euclidean distance from features

  # output name within GRASS GIS
  out_euclidean <- paste0(x, "_inf_nearest")
  # if there is not transformation
  if(is.null(transform)) {
    # if the user provides an output map name, overwrite it
    if(!is.null(output_map_name)) out_euclidean <- output_map_name
    # if no transformation, the output influence is the euclidean distance
    out_influence <- out_euclidean
  }

  # region
  rgrass7::execGRASS("g.region", raster = x, flags = flags_region)
  if(!is.null(extent_x_cut))
    rgrass7::execGRASS("g.region", w = extent_x_cut[1], e = extent_x_cut[2], align = x, flags = flags_region)
  if(!is.null(extent_y_cut))
    rgrass7::execGRASS("g.region", s = extent_y_cut[1], n = extent_y_cut[2], align = x, flags = flags_region)

  # print message
  if(!quiet) rgrass7::execGRASS("g.message", message = "Calculating Euclidean distance...")
  # calculate
  rgrass7::execGRASS("r.grow.distance", input = x, distance = out_euclidean, metric = metric, flags = flags)
  # check if it should be deleted afterwards
  # check if the layer already exists first; if not, remove it in the end.
  if(remove_intermediate)
    if(!is.null(transform)) to_remove <- c(to_remove, out_euclidean)

  # transform distance
  if(!is.null(transform)) {

    if(transform == "log") {

      # log distance
      # message
      out_message <- "Calculating log-distance..."
      # expression
      expression_influence <- sprintf("log(A + %f, %f)", dist_offset, log_base)
      if(print_expression) print(expression_influence)
    }

    if(transform == "sqrt") {

      # sqrt distance
      # message
      out_message <- "Calculating sqrt-distance..."
      # expression
      expression_influence <- sprintf("sqrt(A + %f)", dist_offset)
      if(print_expression) print(expression_influence)

    }

    if(transform == "exp_decay") {
      # define zoi or half life, depending on which is given as input

      # if zoi is given, it is used.
      # if zoi is not given:
      if(is.null(zoi)) {
        # and half_life is not given either:
        if(is.null(half_life)) {
          # lambda is calculated from exp_decay_parms
          lambda <- exp_decay_parms[2]
        } else {
          # if half_life is given, lambda is calculated from that
          zoi <- half_life * zoi_hl_ratio # not used!
          lambda <- log(2)/half_life
        }
      } else {
        # if zoi is given, it is used to
        # define half_life and lambda
        half_life <- zoi/zoi_hl_ratio
        lambda <- log(2)/half_life
      }

      # exponential decay influence
      # message
      out_message <- "Calculating exponential decay influence..."
      # expression
      expression_influence <- sprintf("%f * exp(-%f * A)", exp_decay_parms[1], lambda)
      if(print_expression) print(expression_influence)

    }

    if(transform == "bartlett") {

      # betlett (tent-shaped or linear decay) influence
      # message
      out_message <- "Calculating Bartlett influence..."
      # expression
      expression_influence <- sprintf("if(A <= %f, 1 - (1/%f) * A, 0)", zoi, zoi)
      if(print_expression) print(expression_influence)

    }

    if(transform == "Gauss" | transform == "half_norm") {
      # define zoi or half life, depending on which is given as input

      # if zoi is given, it is used.
      # if zoi is not given:
      if(is.null(zoi)) {
        # and half_life is not given either:
        if(is.null(half_life)) {
          # lambda is calculated from hnorm_decay_parms
          lambda <- 0.5/(hnorm_decay_parms[2]**2) # lambda = 0.5/(sigma^2)
        } else {
          # if half_life is given, lambda is calculated from that
          zoi <- half_life * sqrt(zoi_hl_ratio) # not used!
          lambda <- log(2)/(half_life**2)
        }
      } else {
        # if zoi is given, it is used to
        # define half_life and lambda
        half_life <- zoi/sqrt(zoi_hl_ratio)
        lambda <- log(2)/(half_life**2)
      }

      # exponential decay influence
      # message
      out_message <- "Calculating half-normal decay influence..."
      # expression
      expression_influence <- sprintf("%f * exp(-%f * pow(A, 2))", hnorm_decay_parms[1], lambda)
      if(print_expression) print(expression_influence)

    }

    if(transform %in% c("threshold", "step")) {

      # threshold influence
      # message
      out_message <- "Calculating threshold influence..."
      # expression
      expression_influence <- sprintf("if(A < %f, %f, 0)", zoi, constant_influence)
      if(print_expression) print(expression_influence)
    }

    # if the user provides an output map name, use it
    if(!is.null(output_map_name))
      out_influence <- output_map_name
    else {
      # otherwise, name it according to the method
      out_influence <- paste0(x, "_inf_nearest_", transform)
      # and maybe the zoi
      zoi_methods <- c("exp_decay", "bartlett", "Gauss", "half_norm", "threshold", "step")
      if(!is.null(transform))
        if(transform %in% zoi_methods) out_influence <- paste0(out_influence, zoi)
    }

    # print message
    if(!quiet) rgrass7::execGRASS("g.message", message = out_message)
    # set region
    rgrass7::execGRASS("g.region", raster = out_euclidean, flags = flags_region)
    # calculate
    rgrass7::execGRASS("r.mapcalc.simple", expression = expression_influence,
                       a = out_euclidean, output = out_influence, flags = flags)
  }

  # return only names
  return(out_influence)
}


# backup
# calc_influence_nearest <- function(
#   x,
#   zoi = NULL,
#   transform = NULL, #c("log", "sqrt", "exp_decay", "bartlett", "Gauss", "half_norm")[1],
#   where = c("R", "GRASS")[1],
#   log_base = exp(1),
#   zoi_hl_ratio = 4,
#   half_life = NULL,
#   exp_decay_parms = c(1, 0.01),
#   dist_offset = 1,
#   extent_x_cut = NULL,
#   extent_y_cut = NULL,
#   plotit = FALSE,
#   output_map_name = NULL,
#   metric = c("euclidean", "geodesic", "squared", "maximum", "manhattan")[1],
#   remove_intermediate = TRUE,
#   print_expression = TRUE,
#   quiet = TRUE, overwrite = FALSE,
#   ...) {
#
#   # Run in R
#   if(where %in% c("R", "r")) {
#     if(is.null(extent_x_cut)) extent_x_cut <- terra::ext(x)[c(1,2)]
#     if(is.null(extent_y_cut)) extent_y_cut <- terra::ext(x)[c(3,4)]
#
#     # check if the input is a terra or raster object
#     if(class(x) %in% c("SpatRaster")) {
#       use_terra <- TRUE
#     } else {
#       if(class(x) %in% c("RasterLayer", "RasterBrick", "RasterStack")) {
#         use_terra <- FALSE
#       } else {
#         classes <- c("SpatRaster", "RasterLayer", "RasterBrick", "RasterStack")
#         stop(paste0("Please make sure x is an object of one of these classes: ",
#                     paste(classes, collapse = ","), "."))
#       }
#     }
#
#     # Euclidean distance
#     dist_r <- terra::distance(x, ...)
#
#     # transform distance
#     if(!is.null(transform))
#       if(transform == "log") dist_r <- log(dist_r+dist_offset, base = log_base) else
#         if(transform == "sqrt") dist_r <- sqrt(dist_r+dist_offset) else
#           if(transform == "exp_decay") {
#             # define zoi or half life, depending on which is given as input
#
#             # if zoi is given, it is used.
#             # if zoi is not given:
#             if(is.null(zoi)) {
#               # and half_life is not given either:
#               if(is.null(half_life)) {
#                 # lambda is calculated from exp_decay_parms
#                 lambda <- exp_decay_parms[2]
#               } else {
#                 # if half_life is given, lambda is calculated from that
#                 zoi <- half_life * zoi_hl_ratio # not used!
#                 lambda <- log(2)/half_life
#               }
#             } else {
#               # if zoi is given, it is used to
#               # define half_life and lambda
#               half_life <- zoi/zoi_hl_ratio
#               lambda <- log(2)/half_life
#             }
#
#             dist_r <- exp_decay_parms[1] * exp(-lambda * dist_r)
#           } else
#             if(transform == "bartlett") {
#               dist_r <- (1 - (1/zoi)*dist_r)
#               zero <- dist_r
#               size_landscape <- prod(dim(zero)[c(1:2)])
#               values(zero) <- rep(0, size_landscape)
#               dist_zero_stack <- c(dist_r, zero)
#
#               if(use_terra) {
#                 dist_r <- terra::app(dist_zero_stack, "max")
#               } else {
#                 dist_zero_stack <- raster::stack(dist_zero_stack)
#                 dist_r <- raster::calc(dist_zero_stack, max)
#               }
#
#             } else
#               if(transform %in% c("Gauss", "half_norm")) {
#                 # define zoi or half life, depending on which is given as input
#
#                 # if zoi is given, it is used.
#                 # if zoi is not given:
#                 if(is.null(zoi)) {
#                   # and half_life is not given either:
#                   if(is.null(half_life)) {
#                     # lambda is calculated from exp_decay_parms
#                     lambda <- 0.5/(hnorm_decay_parms[2]**2) # lambda = 0.5/(sigma^2)
#                   } else {
#                     # if half_life is given, lambda is calculated from that
#                     zoi <- half_life * sqrt(zoi_hl_ratio) # not used!
#                     lambda <- log(2)/(half_life**2)
#                   }
#                 } else {
#                   # if zoi is given, it is used to
#                   # define half_life and lambda
#                   half_life <- zoi/sqrt(zoi_hl_ratio)
#                   lambda <- log(2)/(half_life**2)
#                 }
#
#                 dist_r <- hnorm_decay_parms[1] * exp(-lambda * (dist_r**2))
#               } else
#                 stop("You should select an appropriate transformation method for distance.")
#
#     # rename nearest influence layer
#     names(dist_r) <- "influence_nearest"
#     # should the result be plotted?
#     if(plotit) plot(dist_r)
#
#     # return cropped distance layer
#     if(use_terra)
#       inf_nearest <- terra::crop(dist_r, terra::ext(c(extent_x_cut, extent_y_cut)))
#     else
#       inf_nearest <- raster::crop(dist_r, raster::extent(c(extent_x_cut, extent_y_cut)))
#
#     return(inf_nearest)
#   } else {
#
#     # Run in GRASS GIS
#     if(where %in% c("GRASS", "grass", "GRASS GIS", "grass gis")) {
#       inf_nearest <- calc_influence_nearest_GRASS(x = x,
#                                                   zoi = zoi,
#                                                   transform = transform,
#                                                   log_base = log_base,
#                                                   zoi_hl_ratio = zoi_hl_ratio,
#                                                   half_life = half_life,
#                                                   exp_decay_parms = exp_decay_parms,
#                                                   hnorm_decay_parms = hnorm_decay_parms,
#                                                   dist_offset = dist_offset,
#                                                   extent_x_cut = extent_x_cut,
#                                                   extent_y_cut = extent_y_cut,
#                                                   plotit = plotit,
#                                                   output_map_name = output_map_name,
#                                                   metric = metric,
#                                                   remove_intermediate = remove_intermediate,
#                                                   print_expression = print_expression,
#                                                   quiet = quiet,
#                                                   overwrite = overwrite,
#                                                   ...)
#
#       return(inf_nearest)
#     }
#   }
#
# }
