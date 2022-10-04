#' Calculates the zone of influence from the nearest feature
#' and the cumulative zone of influence of multiple features
#'
#' This function takes in a raster with locations of infrastructure and calculates
#' either (1) a raster representing the Zone of Influence (ZoI) of the neareast feature or (2)
#' a raster representing the cumulative Zone of Influence of multiple features, or both.
#' This function takes in a raster with locations or counts of
#' infrastructure and calculates a raster (or set of rasters, in case there is
#' more the one value for `radius`) representing either of two zone of influence (ZoI)
#' metrics for type of infrastructure: (1) the ZoI of the neareast feature, (2) the cumulative
#' ZoI of multiple features, or (3) both ZoI metrics. Zones of influence
#' are defined by functions that decay with the distance from each
#' infrastructure and their rate of decay is controlled by the ZoI radius
#' (`radius`), which defines how far the influence of an infrastructure
#' feature goes. To see more information on each ZoI metric, see
#' [oneimpact::calc_zoi_nearest()] and [oneimpact::calc_zoi_cumulative()].
#'
#' @param x `[RasterLayer,SpatRaster]` \cr Raster representing locations of features,
#' preferentially with positive value where the features
#' are located and either 0 or NA elsewhere.
#' Alternatively, `x` might be a binary (dummy)
#' spatial variable representing the presence of linear or area features, with
#' 0 or NA/no-data as background.
#' `x` can be a `RasterLayer` from [raster] package or a [SpatRaster] from
#' [terra] package. If `where = "GRASS"`, `x` must be a string corresponding
#' to the name of the input map within a GRASS GIS location and mapset.
#'
#' The default parameters assume that the input `x` presents zeros as the background
#' value, where infrastructure or disturbance are absent. Therefore, to deal correctly
#' with the computation of both ZoI metrics, by default we set the parameter `zeroAsNA = TRUE`.
#' If, in contrast, the input map `x` has `NA` as background values, the parameter
#' `zeroAsNA` should be set to `FALSE`.
#'
#' @param radius `[numeric(1)]` \cr Radius of the zone of influence (ZoI),
#' the distance at which the ZoI vanishes or goes below a given minimum limit value
#' `zoi_limit`. See [oneimpact::zoi_functions()] for details.
#' It can be a single value or a vector of values, in which case
#' several ZoI layers (one for each radius) are created.
#'
#' @param type `[character(1)="circle"]{"circle", "Gauss", "rectangle",
#' "exp_decay", "bartlett", "threshold", "step"}` \cr
#' Shape of the zone of influence. See [oneimpact::calc_zoi_nearest()] for
#' details.
#'
#' @param zoi_metric `[character(1)="all"]{"all", "nearest", "cumulative"}` \cr
#' Which metric of zone of influence should be computed. Either `"all"`, `"nearest"`,
#' or `"cumulative"`.
#'
#' @param output_type `[character(1)="cumulative_zoi"]{"cumulative_zoi",
#' "density"}` \cr
#' For the cumulative ZoI, if `output_type = "cumulative_zoi"` (default), the ZoI weight
#' matrix not not normalized, i.e. the maximum value of the weight matrix at the
#' central pixel value is always 1. This means the values of the input map are
#' summed (considering a decay with distance within the neighborhood) and the
#' output map presents values higher than 1. If `output_type = "density"`, the
#' weight matrix is normalized before the filtering process, leading to values
#' in the outmap map generally lower than 1. This parameter is ignored for
#' the ZoI of the nearest feature.
#'
#' @returns If the calculations are performed in R (`where = "R"`),
#' a `RasterLayer`/`RasterStack` or [SpatRaster] object
#' (according to the input `x` map)
#' with the either the zone of influence of the nearest feature
#' (if `zoi_metric = "nearest"`), the cumulative zone of influence of multiple
#' features (if `zoi_metric = "cumulative"`), or both metrics
#' (if `zoi_metric = "all"`, the default). \cr
#' If the computation is done in GRASS GIS, the output is the name of
#' the output raster map(s) within the GRASS GIS location and mapset of the
#' current session. The user can retrieve these maps to R using
#' [rgrass7::read_RAST()] or export them outside GRASS using the
#' `r.out.gdal` module, for instance.
#'
#' @seealso Fore more details on each of the ZoI metrics,
#' other function parameters, and their specific details, see
#' [oneimpact::calc_zoi_nearest()] and [oneimpact::calc_zoi_cumulative()].
#'
#' @example examples/calc_zoi_example.R
#'
#' @export
calc_zoi <- function(x,
                     radius = radius,
                     type = type,
                     zoi_metric = c("all", "nearest", "cumulative")[1],
                     where = c("R", "GRASS")[1],
                     zeroAsNA = TRUE,
                     output_type = c("cumulative_zoi", "density")[1],
                     ...) {

  # check if the input is a terra or raster object
  if(where == "R")
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

  # nearest influence
  if(zoi_metric %in% c("all", "near", "nearest", "dist", "distance")) {
    nearest_r <- calc_zoi_nearest(x = x,
                                  radius = radius,
                                  type = type,
                                  zeroAsNA = zeroAsNA,
                                  ...)
  }

  # cumulative influence
  if(zoi_metric %in% c("all", "cum", "cumulative", "density")) {
    cumulative_r <- calc_zoi_cumulative(x = x,
                                        radius = radius,
                                        type = type,
                                        zeroAsNA = !zeroAsNA,
                                        output_type = output_type,
                                        ...)
  }

  # stack
  if(zoi_metric == "all") {
    if(where == "R") {
      if(use_terra) r_stk <- do.call(c, list(nearest_r, cumulative_r)) else
        r_stk <- raster::stack(nearest_r, cumulative_r)
    } else {
      r_stk <- c(nearest_r, cumulative_r)
    }
  } else {
    if(zoi_metric %in% c("near", "nearest", "dist", "distance")) {
      r_stk <- nearest_r
    } else {
      r_stk <- cumulative_r
    }
  }

  r_stk
}
