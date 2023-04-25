#' Calculate the cumulative zone of influence of multiple features
#'
#' @description This function takes in a raster with locations or counts of
#' infrastructure
#' and calculates a raster (or set of rasters, in case there is more the one
#' value for `radius`) representing the cumulative zone of influence (ZOI)
#' or density of features in space. The process is done through a spatial
#' filter/moving window/neighborhood analysis. The zones of influence
#' (or weight matrices) are defined by functions that decay with the distance
#' from each infrastructure feature and their rate of decay is controlled by the
#' ZOI radius (`radius`), which defines how far the influence of an infrastructure
#' feature goes. By using moving window analyses, the cumulative ZOI account for
#' the sum of the influence of multiple features across space, weighted by
#' the distance to these features according the the ZOI shape.
#' For more details on the ZOI functions, see [oneimpact::zoi_functions()].
#'
#' The procedure might be computed in both R and GRASS GIS. In R, the
#' neighborhood analysis is done with the [terra::focal()] function. In GRASS,
#' different modules might be used for the computation: `r.resamp.filter`,
#' `r.mfilter`, or `r.neighbors`. See details for their differences. In GRASS, it
#' requires an active connection between the R session and a GRASS GIS
#' location and mapset (through the package [rgrass]), and that the input
#' maps are already loaded within this GRASS GIS mapset.
#' If the calculations are done in R, the input is a (set of) raster map(s)
#' and the function returns another (set of) raster map(s). If the calculations
#' are done within GRASS GIS, the input is the name of a raster map already
#' loaded in a GRASS GIS location and mapset, and the function returns
#' only the name of the output map. This map is stored in the the GRASS GIS
#' location/mapset, and might be retrieved to R through the
#' [rgrass::read_RAST()] function or exported outside GRASS using the
#' `r.out.gdal` module, for instance.
#'
#' @details
#' The input raster is supposed to
#' represent the location of point, line, or polygon infrastructure
#' (e.g. houses, roads, mining areas), but any landscape variable whose
#' representation might be one of those would fit here
#' (e.g. areas of forest or any other habitat type or land cover).
#' We recommend that the input raster has a metric projection, so that distances
#' and zones of influence are based on distance to infrastructure measured in meters.
#'
#' ## Zone of Influence functions and weight matrices
#'
#' The neighborhood analysis to define the cumulative ZOI can be
#' computed with different functions/filters. The options currently implemented
#' are:
#' - circular/threshold matrix: the circular filter (`type = "circle"` or
#' `type = "threshold"` or `type = "step"`) is a matrix with constant weights
#' in which the parameter `radius` corresponds to the radius of the circle
#' centered on the central pixel. It is similar to a circular buffer matrix.
#' - Gaussian matrix: the Gaussian filter (`type = "Gauss"` or `type = "gauss"`
#' or `type = "gaussian_decay"`) is a matrix with weights following a Gaussian
#' or Normal decay. The Gaussian curve is 1 at the central cell and
#' is parameterized on the `radius` and
#' the `zoi_limit`, which controls how fast the curve decreases with distance.
#' See [oneimpact::zoi_functions()] for details.
#' - Exponential decay matrix: the exponential decay filter
#' (`type = "exp_decay"`) is a matrix with weights following an exponential
#' decay curve, with value 1 in the central cell and
#' parameterized on the `radius` and the `zoi_limit`.
#' See [oneimpact::zoi_functions()] for details.
#' - Rectangular matrix: the rectangular filter (`type = "rectangle"`
#' or `type = "box"`) is
#' a weight matrix whose shape is a square of dimensions \eqn{n} x \eqn{n},
#' with \eqn{n = 2 * radius}.
#' - Bartlett or linear decay matrix: the Bartlett, linear, or tent decay filter
#' (`type = "bartlett"` or `type = "linear_decay"` or `type = "tent_decay"`)
#' is a weight matrix whose value is 1 in the central cell and whose weights
#' decrease linearly up to zero at a distance equals `radius`.
#' See [oneimpact::zoi_functions()] for details.
#' - user-customized filter: if `type = "mfilter"`, `radius` is not
#' numeric but should be a user-defined matrix of weights. Examples are ones
#' created through [oneimpact::filter_create()], [terra::focalMat()],
#' [smoothie::kernel2dmeitsjer()], or matrices created by hand.
#'
#' Weight matrices might differ from the expected decay function depending on
#' the intended resolution - the finer the resolution, the more detailed and
#' correspondent to the original functions the matrix will be.
#'
#' ## Algorithms in GRASS GIS
#'
#' In GRASS GIS, different modules might be used for the computation,
#' `r.resamp.filter`, `r.mfilter`, or `r.neighbors`. The module to be used is
#' controlled by the parameter `g_module`. These algorithms provide different
#' capabilities and flexibility.
#' - `r.resamp.filter` seems to be the fastest one
#' in most cases, but has less flexibility in the choice of the zone of influence
#' function. The algorithm calculates the weighted density of features, which
#' might be rescaled to the cumulative ZOI if the appropriate scaling factor
#' (calculated from the weight matrix) is provided. Currently only the
#' filters `type = "bartlett"` and `type = "box"` are implemented. More
#' information about the algorithm
#' [here](https://grass.osgeo.org/grass80/manuals/r.resamp.filter.html).
#' - `r.mfilter` is slower than `r.resamp.filter` but much faster than
#' `r.neighbors`, and allows a flexible choice of the shape of the zone of
#' influence (the wight matrix shape). `r.mfilter` is then the most indicated
#' in terms of a balance between flexibility in the choice of the ZOI shape
#' and computation efficiency.
#' The only inconvenient of `r.mfilter` is that it
#' creates an edge effect with no information in the outer cells of a raster
#' (the number of cells correspond to `radius` or half the size of the weight
#' matrix), so if it is used the users should add a buffer area
#' >= `radius` around the input raster map, to avoid such edge effects.
#' See \url{https://github.com/OSGeo/grass/issues/2184} for more details.
#' - `r.neighbors` is considerably slower than the other algorithms (from 10 to
#' 100 times), but allows a flexible choice of the ZOI shape. Contrary to
#' `r.resamp.filter` and `r.mfilter`, which can only perform a sum of pixel
#' values weighted by the input filter or ZOI, `r.neighbors` might
#' calculate many other statistical summaries within the window of analysis,
#' such as mean, median, standard deviation etc.
#'
#' @param x `[RasterLayer,SpatRaster,character]` \cr Raster representing
#' locations of features, preferentially with positive values where the features
#' are located (or counts of features within a pixel) and 0 elsewhere.
#' Alternatively, `x` might be a binary (dummy) spatial variable representing
#' the presence of linear or area features, with zero as background.
#' `x` can be a `RasterLayer` from [raster] package or a [SpatRaster] from
#' [terra] package. If `where = "GRASS"`, `x` must be a string corresponding
#' to the name of the input map within a GRASS GIS location and mapset.
#' Continuous or discrete raster maps with multiple categories can be binarized
#' to be used as input for `calc_zoi_cumulative()` through
#' [landscapetools::util_binarize()] in R or [oneimpact::grass_binarize()]
#' in GRASS GIS, or through common raster algebra in both
#' environments.
#'
#' Notice that, different from [oneimpact::calc_zoi_nearest()], the input maps `x`
#' must have zero as background values, instead of NA. In R it is possible to
#' account for NA background values by setting `zeroAsNA = TRUE` for the computation
#' of the cumulative ZOI.
#' In GRASS, maps without NA as background might be changed into maps with 0 as
#' background to be used in  `calc_zoi_cumulative()`
#' through [raster algebra](https://grass.osgeo.org/grass78/manuals/r.mapcalc.html)
#' and e.g. through the use of the module
#' [`r.null`](https://grass.osgeo.org/grass80/manuals/r.null.html).
#'
#' @param radius `[numeric(1)=100]` \cr Radius or scale of the moving
#' window for neighborhood analysis, used to calculate the cumulative zoi and
#' density. The radius represent
#' the distance at which the ZOI vanishes or goes below a given minimum limit value
#' `zoi_limit`. It can be a single value or a vector of values, in which case
#' several cumulative ZOI or density maps (one for each radius) are created.
#' For `type = "circle"`, the `radius` corresponds to the radius of the
#' circle filter. For `type = "Gauss"` and `type = "exp_decay"`, `radius`
#' corresponds to the distance where the Gaussian or exponential decay function
#' decrease or a small `zoi_limit` value. If `type = "bartlett"`, `radius`
#' is the distance at which the filter reaches zero, after a linear decay
#' from the central pixel. If `type = "rectangle"`, `radius`
#' corresponds to half the size of the side of a square filter.
#' If `type = "mfilter"`, radius is not a numeric value but a matrix itself,
#' defined by the user. See description in the details and [oneimpact::zoi_functions()]
#' for more details on the zone of influence function shapes and parameters.
#'
#' @param type `[character(1)="circle"]{"circle", "Gauss", "rectangle",
#' "exp_decay", "bartlett", "threshold", "step", "mfilter"}` \cr
#' Type of filter (shape of the zone of influence) used to calculate the
#' cumulative ZOI or density. See details.
#'
#' @param zoi_limit `[numeric(1)=0.05]` \cr For non-vanishing functions
#' (e.g. `exp_decay`, `gaussian_decay`), this value is used to set the relationship
#' between the ZOI radius and the decay functions:
#' `radius` is defined as the minimum distance at which the ZOI assumes values
#' below `zoi_limit`. The default is 0.05. This parameter is used only
#' if `radius` is not `NULL`.
#'
#' @param output_type `[character(1)="cumulative_zoi"]{"cumulative_zoi",
#' "density"}` \cr If `output_type = "cumulative_zoi"` (default), the ZOI weight
#' matrix not not normalized, i.e. the maximum value of the weight matrix at the
#' central pixel value is always 1. This means the values of the input map are
#' summed (considering a decay with distance within the neighborhood) and the
#' output map presents values higher than 1. If `output_type = "density"`, the
#' weight matrix is normalized before the filtering process, leading to values
#' in the outmap map generally lower than 1.
#'
#' @param where `[character(1)="R"]{"R", "GRASS"}` \cr Where should the
#' computation be done? Default is `"R"`. If `where = "GRASS"`, the R session
#' must be linked to an open GRASS GIS session in a specific location and mapset.
#'
#' @param g_module `[character(1)="r.mfilter"]{"r.mfilter",
#' "r.resamp.filter", "r.neighbors"}` \cr
#' If `where = "GRASS"`, which algorithm should be used to compute the cumulative
#' ZOI? See details for their description.
#'
#' @param min_intensity `[numeric(1)=0.01]` \cr Minimum intensity of the
#' exponential and Gaussian decay functions to
#' define the radius of the window that define the filter. See
#' [oneimpact::filter_create()] for details.
#' @param max_dist `[numeric(1)=50000]` \cr Maximum size (in meters) to
#' define the radius of the window that defines the filter. Only
#' applicable for exponential and Gaussian decay functions. See
#' [oneimpact::filter_create()] for details.
#'
#' @param zeroAsNA `[logical(1)=FALSE]` \cr If `TRUE` treats cells that are
#' `NA` as if they were zero.
#'
#' @param extent_x_cut,entent_y_cut `[numeric vector(2)=c(0,1)]` \cr Vector
#' representing the minimum and
#' maximum extent in x and y for the final output, in the format c(min,max).
#' It is intended to keep only a region of interest but consider the
#' surroundings when calculating the cumulative ZOI or density. This might be
#' especially useful for example in the use of the `r.mfilter` algorithm in
#' GRASS, in which the edges of the region are excluded from the computation.
#' The default is to keep the same extent of the input raster.
#' @param na.policy `[character(1)="omit"]` \cr Can be used to determine the
#' cells of `x` for which focal values should be computed. Must be one of "all"
#' (compute for all cells), "only" (only for cells that are NA) or "omit"
#' (skip cells that are NA). Note that the value of this argument does not
#' affect which cells around each focal cell are included in the computations
#' (use na.rm=TRUE to ignore cells that are NA for that). See [terra::focal()]
#' for details. Only used when `where = "R"`.
#' @param na.rm `[logical(1)=FALSE]` \cr Should missing values be removed for
#' filtering calculations? Option for the neighborhood analysis performed
#' through the [terra::focal()] function. Only used when `where = "R"`.
#' @param ... Other arguments to be used within [oneimpact::filter_create()]
#' or [terra::focal()].
#'
#' @param g_output_map_name `[character(1)=NULL]` \cr Name of the output map. Only
#' used when `where = "GRASS"`. If `NULL` (default), a standard name is created
#' based on the name of the input map `x`, the ZOI shape `type`, and the ZOI
#' radius `radius`.
#' @param g_parallel `[logical(1)=TRUE]` \cr If the `r.mfilter` tool is used for
#' computing the cumulative ZOI in GRASS GIS, `g_parallel` tells if parallelization
#' should be used. Default is `TRUE`. The option is ignored for other GRASS tools
#' (parameter `g_module`) and for computation in R (`where = "R"`).
#' @param g_input_as_region `[logical(1)=TRUE]` \cr Should the input map `x` be
#' used to redefine the working GRASS region before cumulative ZOI calculation?
#' If `TRUE`, `x` is used to define the region with `g.region`. If `FALSE`,
#' the region previously defined in the GRASS GIS session is used for computation.
#' @param g_remove_intermediate `[logical(1)=TRUE]` \cr Should the intermediate
#' maps created for computing the output map be excluded in the end of the
#' process? Only used when `where = "GRASS"`.
#' @param g_overwrite `[logical(1)=FALSE]` \cr If the a map already exists with the
#' name `g_output_map_name` in the working GRASS GIS location and mapset, should
#' it be overwritten? Only used when `where = "GRASS"`.
#' @param verbose `[logical(1)=FALSE]` \cr Should messages of the computation steps
#' be printed in the prompt along the computation?
#'
#' @returns If the calculations are performed in R (`where = "R"`),
#' a `RasterLayer` or [SpatRaster] (according to the input `x` map)
#' with the cumulative zone of influence or density of features. While the
#' cumulative ZOI uses a ZOI/weight matrix rescaled to 1 at the central pixel
#' (creating values in the output map which might go well beyond 1), the
#' density of features uses a normalized ZOI/weight matrix (with all values
#' summing 1), what created values smaller than one in the output map.
#' if multiple `radius` values are given, a `RasterBrick` or multi-layer
#' `SpatRaster`, with the cumulative ZOI or density maps for each ZOI radius. \cr
#' If the computation is done in GRASS GIS, the output is the name of
#' the output raster map(s) within the GRASS GIS location and mapset of the
#' current session. The user can retrieve these maps to R using
#' [rgrass::read_RAST()] or export them outside GRASS using the
#' `r.out.gdal` module, for instance.
#'
#' @seealso See [oneimpact::zoi_functions()] for some ZOI function shapes and
#' [oneimpact::filter_create()] for options to create weight matrices. \cr
#' See also [smoothie::kernel2dmeitsjer()], [terra::focalMat()], and
#' [raster::focalWeight()] for other functions to create filters or weight matrices. \cr
#' See
#' [r.mfilter](https://grass.osgeo.org/grass80/manuals/r.mfilter.html),
#' [r.resamp.filter](https://grass.osgeo.org/grass80/manuals/r.resamp.filter.html), and
#' [r.neighbors](https://grass.osgeo.org/grass80/manuals/r.neighbors.html) for
#' GRASS GIS implementations of neighborhood analysis.\cr
#' See [oneimpact::calc_zoi_nearest()] for the computation of the zone of influence
#' of the nearest feature only.
#'
#' @example examples/calc_zoi_cumulative_example.R
#' @example examples/calc_zoi_cumulative_grass_example.R
#'
#' @export
calc_zoi_cumulative <- function(
  x,
  radius = 100,
  type = c("circle", "Gauss", "rectangle", "exp_decay",
           "bartlett", "threshold", "mfilter")[1],
  where = c("R", "GRASS")[1],
  output_type = c("cumulative_zoi", "density")[1],
  zoi_limit = 0.05,
  min_intensity = 0.01,
  max_dist = 50000,
  zeroAsNA = FALSE,
  extent_x_cut = NULL,
  extent_y_cut = NULL,
  na.policy = "omit",
  na.rm = TRUE,
  g_module = c("r.mfilter", "r.resamp.filter", "r.neighbors")[1],
  g_parallel = TRUE,
  g_output_map_name = NULL,
  g_input_as_region = FALSE,
  g_remove_intermediate = TRUE,
  g_overwrite = FALSE,
  verbose = FALSE, ...) {

  # Run in R
  if(where %in% c("R", "r")) {
    if(is.null(extent_x_cut)) extent_x_cut <- terra::ext(x)[c(1,2)]
    if(is.null(extent_y_cut)) extent_y_cut <- terra::ext(x)[c(3,4)]

    zoi_cumulative <- calc_zoi_cumulative_r(x,
                                            radius = radius,
                                            type = type,
                                            output_type = output_type,
                                            zoi_limit = zoi_limit,
                                            min_intensity = min_intensity,
                                            max_dist = max_dist,
                                            zeroAsNA = zeroAsNA,
                                            extent_x_cut = extent_x_cut,
                                            extent_y_cut = extent_y_cut,
                                            na.policy = na.policy,
                                            na.rm = na.rm,
                                            verbose = verbose,
                                            ...)

    return(zoi_cumulative)
  } else {

    # Run in GRASS GIS
    if(where %in% c("GRASS", "grass", "GRASS GIS", "grass gis")) {
      zoi_cumulative <- calc_zoi_cumulative_grass(x = x,
                                                  radius = radius,
                                                  type = type,
                                                  output_type = output_type,
                                                  g_module = g_module,
                                                  zoi_limit = zoi_limit,
                                                  min_intensity = min_intensity,
                                                  max_dist = max_dist,
                                                  extent_x_cut = extent_x_cut,
                                                  extent_y_cut = extent_y_cut,
                                                  zeroAsNA = zeroAsNA,
                                                  g_parallel = g_parallel,
                                                  g_output_map_name = g_output_map_name,
                                                  g_input_as_region = g_input_as_region,
                                                  g_remove_intermediate = g_remove_intermediate,
                                                  g_overwrite = g_overwrite,
                                                  verbose = verbose,
                                                  ...)

      return(zoi_cumulative)
    }
  }


}

# implementation in R
calc_zoi_cumulative_r <- function(
  x,
  radius = 100,
  type = c("circle", "Gauss", "rectangle", "exp_decay",
           "bartlett", "threshold", "mfilter")[1],
  output_type = c("cumulative_zoi", "density")[1],
  zoi_limit = 0.05,
  zoi_hl_ratio = NULL,
  half_life = NULL,
  exp_decay_parms = c(1, 0.01),
  min_intensity = 0.01,
  max_dist = 50000,
  zeroAsNA = FALSE,
  extent_x_cut = terra::ext(x)[c(1,2)],
  extent_y_cut = terra::ext(x)[c(3,4)],
  na.policy = "omit",
  na.rm = TRUE,
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

  # check if the input raster presents only a single value (1,NA)
  # if so, transform it into a binary map (1,0)
  r0 <- x
  if(zeroAsNA) {
    if(use_terra) {
      r0 <- terra::classify(r0, cbind(NA, 0)) # binary map
    } else {
      r0 <- raster::reclassify(r0, cbind(NA, 0)) # binary map
    }
  }
  # plot(r0)

  # check filter names
  if(type %in% c("gauss", "gaussian_decay", "normal_decay")) type <- "Gauss"
  if(type %in% c("bartlett_decay", "linear_decay", "tent_decay")) type <- "bartlett"
  if(type %in% c("threshold", "threshold_decay", "step_decay")) type <- "circle"

  # cum zoi vs density
  if(output_type %in% c("cumulative_zoi", "zoi", "cumulative")) normalize <- FALSE
  else if(output_type %in% c("density", "Density")) normalize <- TRUE

  # define filters
  if(type %in% c("exp_decay", "bartlett", "circle", "threshold", "rectangle", "Gauss")) {
    if(length(radius) == 1) {
      filt <- oneimpact::filter_create(r0, radius = radius,
                                       type = type, zoi_limit = zoi_limit,
                                       zoi_hl_ratio = zoi_hl_ratio,
                                       half_life = half_life,
                                       max_dist = max_dist,
                                       min_intensity = min_intensity,
                                       normalize = normalize, ...)
    } else {
      filt <- purrr::map(radius, function(z, ...) {
        oneimpact::filter_create(r0, radius = z, type = type,
                                 zoi_limit = zoi_limit,
                                 zoi_hl_ratio = zoi_hl_ratio,
                                 half_life = half_life,
                                 max_dist = max_dist,
                                 min_intensity = min_intensity,
                                 normalize = normalize, ...)
      })
    }
  }

  if(type == "mfilter") {
    filt <- radius
  }

  # those methods were put above with the filter_create function
  # type = c("circle", "rectangle", "threshold", "step")
  # only Gauss is kept here, so far
  # if(type %in% c("circle", "Gauss", "rectangle", "threshold", "step")) {
  # if(type %in% c("Gauss")) {
  #
  #   # if(type %in% c("threshold", "step")) type <- "circle"
  #   # if(type == "rectangle") radius = 2*radius # for this case d is the side of the square
  #
  #   if(length(radius) == 1) {
  #     filt <- terra::focalMat(r0, d = radius, type = type)
  #     if(!normalize) filt <- filt/max(filt, na.rm = T)
  #   } else {
  #     filt <- purrr::map(radius, function(z) {
  #       ft <- terra::focalMat(r0, d = z, type = type)
  #       if(!normalize) ft <- ft/max(ft, na.rm = T)
  #       ft
  #     })
  #   }
  # }

  # neighborhood analysis
  if(type == "mfilter") {
    # more than one matrix
    if("list" %in% class(filt)) {
      cuminf <- purrr::map2(filt, 1:length(radius), function(f, z) {
        if(verbose) print(paste0("Calculating for ZOI n. ", z, "..."))
        terra::focal(r0, w = f, na.policy = na.policy, na.rm = na.rm, ...)
      })
      if(use_terra) cumulative_r <- do.call(c, cuminf) else
        cumulative_r <- raster::stack(cuminf)
    } else {
      #only one matrix
      cumulative_r <- terra::focal(r0, w = filt, na.policy = na.policy, na.rm = na.rm, ...)
    }
  } else {
    if(length(radius) == 1) {
      cumulative_r <- terra::focal(r0, w = filt, na.policy = na.policy, na.rm = na.rm, ...)
    } else {
      cuminf <- purrr::map2(filt, radius, function(f, z) {
        if(verbose) print(paste0("Calculating for ZOI = ", z, "..."))
        terra::focal(r0, w = f, na.policy = na.policy, na.rm = na.rm, ...)
      })
      if(use_terra) cumulative_r <- do.call(c, cuminf) else
        cumulative_r <- raster::stack(cuminf)
    }
  }

  # rename cumulative zoi layer
  if(type == "mfilter") {
    if(!is.list(radius)) name <- "zoi_cumulative" else
      name <- paste0("zoi_cumulative", 1:length(radius))
  } else {
    name <- paste0("zoi_cumulative_", type, radius)
  }

  names(cumulative_r) <- name

  # return cropped raster
  if(use_terra)
    terra::crop(cumulative_r, terra::ext(c(extent_x_cut, extent_y_cut)))
  else
    raster::crop(cumulative_r, raster::extent(c(extent_x_cut, extent_y_cut)))
}

# implementation in GRASS
calc_zoi_cumulative_grass <- function(
  x,
  radius = 100,
  type = c("circle", "Gauss", "rectangle", "exp_decay", "bartlett", "threshold", "step", "mfilter")[1],
  g_module = c("r.mfilter", "r.resamp.filter", "r.neighbors")[1],
  output_type = c("cumulative_zoi", "density")[1],
  zoi_limit = 0.05,
  zoi_hl_ratio = NULL,
  half_life = NULL,
  exp_decay_parms = c(1, 0.01),
  hnorm_decay_parms = c(1, 20),
  min_intensity = 0.01,
  max_dist = 50000,
  divisor = 1,
  normalize = FALSE,
  extent_x_cut = NULL,
  extent_y_cut = NULL,
  zeroAsNA = FALSE,
  g_parallel = TRUE,
  g_output_map_name = NULL,
  g_input_as_region = FALSE,
  g_remove_intermediate = TRUE,
  g_overwrite = FALSE,
  verbose = FALSE,
  ...) {

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
  ##### 3. CUT extent to be implemented

  # start by setting the region
  if(g_input_as_region)
    rgrass::execGRASS("g.region", raster = x, flags = flags_region)

  # Differently from the R version, here we do not check if the input map
  # is binary or if it has NA values
  # This should be checked by the users before calculation

  if(zeroAsNA) {
    # print message
    if(verbose) rgrass::execGRASS("g.message", message = "Setting zeros as NULL data...")
    # copy map and use r.null
    set_null_map <- "temp_set_null_input"
    # copy map
    rgrass::execGRASS("g.copy", raster = paste0(x, ",", set_null_map), flags = flags)
    # set nulls
    rgrass::execGRASS("r.null", map = set_null_map, null = "0")
    # define input
    input_bin <- set_null_map
    # check if it should be deleted afterwards
    if(g_remove_intermediate)
      to_remove <- c(to_remove, set_null_map)
  } else {
    input_bin <- x
  }

  # define the name of the output map
  if(!is.null(g_output_map_name)) {
    # given name if this is given as a parameter
    out_map <- g_output_map_name
  } else {
    # define name as input + cumulative + type
    out_map = ifelse(output_type %in% c("cumulative_zoi", "zoi", "cumulative"),
                     paste0(x, "_zoi_cumulative_", type),
                     paste0(x, "_density_", type))
  }

  # Define parameters based on output_type
  # Cumulative ZOI or Density
  if(output_type %in% c("cumulative_zoi", "zoi", "cumulative")) {
    normalize <- FALSE
    message_name <- "cumulative ZOI"
  } else {
    if(output_type %in% c("density", "Density")) {
      normalize <- TRUE
      message_name <- "density"
    }
  }

  # get resolution
  region <- rgrass::gmeta()
  resolution <- region$nsres

  allowed_modules <- c("r.resamp.filter", "r.mfilter", "r.neighbors")
  if(!(g_module %in% allowed_modules))
    stop(paste0("You should use one of the following GRASS GIS modules: ",
                paste(allowed_modules, collapse = ","), "."))

  # perform calculations for "r.resamp.filter"
  if(g_module == "r.resamp.filter") {

    # allowed methods
    resamp_filter_types <- c("box", "bartlett", "gauss", "normal", "hermite",
                             "sinc", "lanczos1", "lanczos2", "lanczos3",
                             "hann", "hamming", "blackman")

    # check methods
    types <- strsplit(type, split = ",")[[1]]

    # check filter names
    for(i in types) {
      if(types[i] %in% c("Gauss", "gaussian_decay", "normal_decay"))
        types[i] <- "gauss"
      if(types[i] %in% c("bartlett_decay", "linear_decay",
                         "Bartlett", "tent_decay"))
        types[i] <- "bartlett"
      if(types[i] %in% c("rectangle"))
        types[i] <- "box"
    }

    # check if filter types are allowed
    if(!all(types %in% resamp_filter_types))
      stop(paste0("For the use of GRASS GIS module 'r.resamp.filter', please ",
                  "choose among the following filters: ",
                  paste(resamp_filter_types, collapse = ","),
                  ". Currently, the following filters were selected: ",
                  paste(types, collapse = ","), "."))

    # define filters for rescaling if we want the cumulative ZOI
    if(output_type %in% c("cumulative_zoi", "zoi", "cumulative")) {

      # here we must implement the filters exactly as they are used in
      # r.resamp.filter

      if(types[1] %in% c("bartlett", "box", "gauss")) {
        filter_count <- radius
        filter_file <- tempfile(paste0("my_filter", filter_count, "_"))
        if(length(radius) == 1) {
          filt <- oneimpact::filter_create(r = resolution, radius = radius,
                                           type = type, zoi_limit = zoi_limit,
                                           zoi_hl_ratio = zoi_hl_ratio,
                                           half_life = half_life,
                                           max_dist = max_dist,
                                           min_intensity = min_intensity,
                                           normalize = TRUE,
                                           save_txt = FALSE,
                                           save_format = "raw",
                                           save_file = filter_file, ...)
        } else {
          filt <- purrr::map2(radius, filter_file, function(z, file, ...) {
            oneimpact::filter_create(r = resolution, radius = z,
                                     zoi_limit = zoi_limit, type = type,
                                     zoi_hl_ratio = zoi_hl_ratio,
                                     half_life = half_life,
                                     max_dist = max_dist,
                                     min_intensity = min_intensity,
                                     divisor = divisor,
                                     normalize = TRUE,
                                     save_txt = FALSE,
                                     save_format = "raw",
                                     save_file = file, ...)
          })
        }
      } else {
        stop("Currently, the cumulative ZOI using 'r.resamp.filter' is only
              implemented for the filters 'bartlett', 'rectangle/box', and
              'gauss'. Please select one of those or change the GRASS GIS
              module to either 'r.mfilter' or 'r.neighbors'.")
      }
    }

    # set parameters for neighborhood analysis

    # only one matrix
    if(length(radius) == 1) {
      out_names <- out_map
      out_names <- paste0(out_names, radius)
      filters <- paste(types, collapse = ",")
      parms <- list(input = input_bin, output = out_names, filter = filters,
                    radius = radius)
    } else {
      # several matrices or zoi values
      if(length(radius) > 1) {
        parms <- purrr::map(radius, function(x)
          list(input = input_bin, output = paste0(out_map, x),
               filter = paste(types, collapse = ","), radius = x))
        out_names <- purrr::map(parms, ~ .$output) |> unlist()
      }
    }

    # run neighborhood analysis
    if(length(radius) > 1) {

      # loop for matrices or zoi values
      for(i in seq_along(radius)) {
        parm <- parms[[i]]
        z <- radius[[i]]

        # change names for creating density map first, in case the
        # final output is the cumulative ZOI
        if(output_type %in% c("cumulative_zoi", "zoi", "cumulative")) {
          out_final <- parm$output
          parm$output <- paste0(parm$output, "_temp")
          if(g_remove_intermediate) to_remove <- c(to_remove, parm$output)
        }

        # message
        msg <- paste0("Calculating ", message_name, " for ", z, ", shape ",
                      type, "...")
        if(verbose) print(msg)
        # region
        # set region
        if(g_input_as_region)
          rgrass::execGRASS("g.region", raster = parm$input, flags = flags_region)
        # calculate
        rgrass::execGRASS(g_module, parameters = parm, flags = flags)

        # rescale to the cumulative ZOI, if this is the intended final output
        if(output_type %in% c("cumulative_zoi", "zoi", "cumulative")) {

          # maximum value from weight matrix
          max_val <- max(filt[[i]], na.rm = TRUE)
          # expression
          expr <- paste0(out_final, " = ", parm$output, "/", max_val)

          # message
          msg <- paste0("Rescaling from density to cumulative ZOI...")
          if(verbose) print(msg)
          print(expr)

          # calculate ZOI map
          out_final <- rgrass::execGRASS("r.mapcalc", expression = expr,
                                          flags = flags)
        }
      }

    } else {

      # for only one matrix/zoi value
      parm <- parms
      z <- radius

      # change names for creating density map first, in case the
      # final output is the cumulative ZOI
      if(output_type %in% c("cumulative_zoi", "zoi", "cumulative")) {
        out_final <- parm$output
        parm$output <- paste0(parm$output, "_temp")
        if(g_remove_intermediate) to_remove <- c(to_remove, parm$output)
      }

      # message
      msg <- paste0("Calculating ", message_name, " for ", z, ", shape ",
                    type, "...")
      if(verbose) print(msg)
      # region
      # set region
      if(g_input_as_region)
        rgrass::execGRASS("g.region", raster = parm$input, flags = flags_region)
      # calculate
      rgrass::execGRASS(g_module, parameters = parm, flags = flags)

      # rescale to the cumulative ZOI, if this is the intended final output
      if(output_type %in% c("cumulative_zoi", "zoi", "cumulative")) {

        # maximum value from weight matrix
        max_val <- max(filt, na.rm = TRUE)
        # expression
        expr <- paste0(out_final, " = ", parm$output, "/", max_val)

        # message
        msg <- paste0("Rescaling from density to cumulative ZOI...")
        if(verbose) print(msg)
        print(expr)

        # calculate ZOI map
        out_final <- rgrass::execGRASS("r.mapcalc", expression = expr,
                                        flags = flags)
      }
    }
  }

  # for r.neighbors, it always performs the average, it is not possible to sum
  # the matrix is always normalized
  # one must use the argument size in conjunction with the weight matrix
  # we must a specific output from save_mfilter for that, with only the numbers
  # r.neighbors input=private_cabins_sub_bin output=test_neighbors method=count size=21 weight=test_neighbors_exp_filt500.txt --o

  # perform calculations for "r.mfilter"
  if(g_module == "r.mfilter") {

    # define filters
    if(type %in% c("exp_decay", "bartlett", "circle", "threshold", "step", "rectangle", "Gauss")) {
      filter_count <- radius
      filter_file <- tempfile(paste0("my_filter", filter_count, "_"))
      if(length(radius) == 1) {
        filt <- oneimpact::filter_create(r = resolution,
                                         radius = radius,
                                         type = type, zoi_limit = zoi_limit,
                                         zoi_hl_ratio = zoi_hl_ratio,
                                         half_life = half_life,
                                         max_dist = max_dist,
                                         min_intensity = min_intensity,
                                         divisor = divisor,
                                         normalize = normalize, save_txt = TRUE,
                                         save_file = filter_file, ...)
      } else {
        filt <- purrr::map2(radius, filter_file, function(z, file, ...) {
          oneimpact::filter_create(r = resolution, radius = z,
                                   zoi_limit = zoi_limit, type = type,
                                   zoi_hl_ratio = zoi_hl_ratio,
                                   half_life = half_life,
                                   max_dist = max_dist,
                                   min_intensity = min_intensity,
                                   divisor = divisor,
                                   normalize = normalize, save_txt = TRUE,
                                   save_file = file, ...)
        })
      }
    }

    # In these cases the matrix is not defined by filter_create
    if(type %in% c("mfilter")) {

      # Filters pre-defined for "mfilter"
      if(type == "mfilter") {
        # set
        filt <- radius
        # normalize if not already normalized, and if they should be
        if(normalize) {
          # only one matrix
          if(is.matrix(filt)) {
            ss <- sum(filt, na.rm = TRUE)
            if(ss != 1) filt <- filt/ss
          } else {
            # if it is a series of matrices
            if(is.list(filt)) {
              ss <- purrr::map(filt, sum, na.rm = TRUE)
              if(any(ss != 1)) filt <- purrr::map2(filt, ss, ~.x/.y)
            }
          }
        }
      }

      # for one matrix only
      if(is.matrix(filt)) {
        filter_count <- 1
        filter_file <- tempfile(paste0("my_filter_", type, filter_count, "_"))
        # save matrix outside R for use within GRASS GIS
        oneimpact::filter_save(filt, radius = "", type = type,
                               divisor = divisor,
                               save_format = c("GRASS_rmfilter"),
                               save_file = filter_file,
                               parallel = g_parallel,
                               separator = " ")
      } else {
        # for multiple matrices
        if(is.list(filt)) {
          filter_count <- 1:length(filt)
          filter_file <- tempfile(paste0("my_filter_", type, filter_count, "_"))
          # save matrices outside R for use within GRASS GIS
          purrr::map2(filt, filter_file, function(f, file, ...) {
            oneimpact::filter_save(f, radius = "", type = type,
                                   divisor = divisor,
                                   save_format = c("GRASS_rmfilter"),
                                   save_file = file,
                                   parallel = g_parallel,
                                   separator = " ")
          })
        }
      }
    }

    # set parameters for neighborhood analysis

    # only one matrix
    if(is.matrix(filt)) {
      out_names <- out_map
      if(type != "mfilter") out_names <- paste0(out_names, radius)
      parms <- list(input = input_bin, output = out_names, filter = filter_file)
    } else {
      # several matrices or zoi values
      if(is.list(filt)) {
        parms <- purrr::map2(filter_count, filter_file, function(x, y)
          list(input = input_bin, output = paste0(out_map, x), filter = y))
        out_names <- purrr::map(parms, ~ .$output) |> unlist()
      }
    }

    # run neighborhood analysis
    if("list" %in% class(filt)) {

      # loop for matrices or zoi values
      for(i in 1:length(filt)) {
        parm <- parms[[i]]
        z <- filter_count[[i]]

        # message
        msg <- paste0("Calculating ", message_name, " for ", z, ", shape ", type, "...")
        if(verbose) print(msg)
        # region
        # set region
        if(g_input_as_region)
          rgrass::execGRASS("g.region", raster = parm$input, flags = flags_region)
        # calculate
        rgrass::execGRASS(g_module, parameters = parm, flags = flags)
      }

    } else {

      # for only one matrix/zoi value
      parm <- parms
      z <- filter_count

      # message
      msg <- paste0("Calculating ", message_name, " for ", z, ", shape ", type, "...")
      if(verbose) print(msg)
      # region
      # set region
      if(g_input_as_region)
        rgrass::execGRASS("g.region", raster = parm$input, flags = flags_region)
      # calculate
      rgrass::execGRASS(g_module, parameters = parm, flags = flags)

    }
  }

  # perform calculations for "r.mfilter"
  if(g_module == "r.neighbors") {
    stop("Usage through 'r.neighbors' not implemented yet.")
  }

  # remove intermediate maps
  remove_flags = ifelse(verbose, "f", c("f", "quiet"))
  if(g_remove_intermediate)
    if(length(to_remove) > 0)
      rgrass::execGRASS("g.remove", type = "rast", name = to_remove,
                         flags = remove_flags)

  # return only names
  return(out_names)
}
