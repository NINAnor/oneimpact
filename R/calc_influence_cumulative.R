#' Calculate cumulative influence of multiple features
#'
#' r.resamp.filter resamples an input raster, filtering the input with an analytic kernel.
#' Each output cell is typically calculated based upon a small subset of the input cells,
#' not the entire input. r.resamp.filter performs convolution (i.e. a weighted sum is
#' calculated for every raster cell).
#'
#' This function takes in a raster with locations of infrastructure and calculates
#' a raster (or set of rasters, in case there is more the one value for `zoi`)
#' representing the density of features in space (through a spatial filter/neighborhood analysis).
#' The neighborhood analysis is done with the [terra::focal()] function.
#'
#' The neighborhood analysis can be done with different methods. The default is a circular filter
#' (`type = "circle"`), in which case the parameter `zoi` corresponds to the radius of the circle
#' centered on the central pixel. Other possibilities are a Gaussian filter
#' (`type = "Gauss"`), in which case scale corresponds to the sigma parameter of a Gaussian
#' functionl and `type = "rectangle"`, in which case the scale corresponds to the
#' size size of the rectangle. For all these methods, these parameters feed the function
#' [terra::focalMat()] to create the input weight matrix. See [terra::focalMat()] for more
#' details.
#'
#' If one wants to use their own filter or weight matrix, it is possible to define `type = "mfilter"`
#' and provide a matrix to the parameter `zoi` instead, such as one created through the
#' [terra::focalMat()] or the [oneimpact::create_filter()] functions.
#'
#' TO IMPROVE1: implement with `terra`. WE SHOULD DETECT IF THE INPUT IS RASTER OR TERRA
#'
#' TO IMPROVE2: do the same in communication with GRASS GIS.
#'
#' @param x `[RasterLayer,SpatRaster]` \cr Raster representing locations of features, with 1 where the features
#' are located and NA elsewhere. Can be a [RasterLayer] from [raster] package or a [SpatRaster] from
#' [terra] package.
#' @param type `[character(1)="circle"]{"circle", "Gauss", "rectangle", "exp_decay", "bartlett", "mfilter"}` \cr
#' Type of filter used to calculate density. See description for details.
#' @param zoi `[numeric(1)=100]` \cr Scale of the neighborhood analysis, used to calculate densities.
#' It can be a single value of a vector of values, in which case several density maps (for each scale)
#' are created. For `type = "circle"`, scale corresponds to the radius of the circle filter. For `type = "Gauss"`,
#' it corresponds to the standard deviation of the Gaussian distribution. If `type = "rectangle"`, it corresponds
#' to the size of the side of a square filter. See [terra::focalMat()] for more details.
#' If `type = "mfilter"`, scale is not a numeric value but a matrix itself, defined by the user. See description
#' in the details.
#' @param extent_x_cut,entent_y_cut `[numeric vector(2)=c(0,1)]` \cr Vector representing the minimum and
#' maximum extent in x and y for the final output, in the format c(min,max). It is intended to keep only
#' a region of interest but consider the surroundings when calculating densities. The default is to
#' keep the same extent of the input raster.
#' @param na.rm `[logical(1)=FALSE]` \cr Should missing values be removed for filtering calculations?
#' Option for the neighborhood analysis, performed through the [terra::focal()] function.
#' @param plotit `[logical(1)=FALSE]` \cr Should the outputs be plotted along the calculation?
#' @param ... Other arguments to be used within [oneimpact::create_filter()] or [terra::focal()].
#'
#' @returns A [RasterLayer] or [SpatRaster] (according to the input `x` map) with the density of features.
#' Alternatively, a `RasterBrick`
#' or multi-layer `SpatRaster`, if multile `zoi` values are given, with the density of features for all scales selected.
#'
#' @example examples/calc_inf_cumulative_example.R
#'
#' @export
calc_influence_cumulative <- function(x,
                                      zoi = 100,
                                      type = c("circle", "Gauss", "rectangle", "exp_decay", "bartlett", "mfilter")[1],
                                      where = c("R", "GRASS")[1],
                                      module = c("r.mfilter", "r.resamp.filter", "r.neighbors")[1],
                                      zoi_hl_ratio = 4,
                                      half_life = NULL,
                                      exp_decay_parms = c(1, 0.01),
                                      min_intensity = 0.01,
                                      max_dist = 50000,
                                      normalize = FALSE,
                                      divisor = 1,
                                      extent_x_cut = NULL,
                                      extent_y_cut = NULL,
                                      na.policy = "omit",
                                      na.rm = TRUE,
                                      plotit = FALSE,
                                      parallel = TRUE,
                                      output_map_name = NULL,
                                      remove_intermediate = TRUE,
                                      overwrite = FALSE,
                                      quiet = TRUE,...) {

  # Run in R
  if(where %in% c("R", "r")) {
    if(is.null(extent_x_cut)) extent_x_cut <- terra::ext(x)[c(1,2)]
    if(is.null(extent_y_cut)) extent_y_cut <- terra::ext(x)[c(3,4)]

    inf_cumulative <- calc_influence_cumulative_r(x,
                                                  zoi = zoi,
                                                  type = type,
                                                  zoi_hl_ratio = zoi_hl_ratio,
                                                  half_life = half_life,
                                                  exp_decay_parms = exp_decay_parms,
                                                  min_intensity = min_intensity,
                                                  max_dist = max_dist,
                                                  normalize = normalize,
                                                  extent_x_cut = extent_x_cut,
                                                  extent_y_cut = extent_y_cut,
                                                  na.policy = na.policy,
                                                  na.rm = na.rm,
                                                  quiet = quiet,
                                                  plotit = plotit, ...)

    return(inf_cumulative)
  } else {

    # Run in GRASS GIS
    if(where %in% c("GRASS", "grass", "GRASS GIS", "grass gis")) {
      inf_cumulative <- calc_influence_cumulative_GRASS(x = x,
                                                        zoi = zoi,
                                                        type = type,
                                                        module = module,
                                                        zoi_hl_ratio = zoi_hl_ratio,
                                                        half_life = half_life,
                                                        exp_decay_parms = exp_decay_parms,
                                                        min_intensity = min_intensity,
                                                        max_dist = max_dist,
                                                        divisor = divisor,
                                                        normalize = normalize,
                                                        extent_x_cut = extent_x_cut,
                                                        extent_y_cut = extent_y_cut,
                                                        parallel = parallel,
                                                        output_map_name = output_map_name,
                                                        remove_intermediate = remove_intermediate,
                                                        overwrite = overwrite,
                                                        quiet = quiet,
                                                        ...)

      return(inf_cumulative)
    }
  }


}

calc_influence_cumulative_r <- function(
  x,
  zoi = 100,
  type = c("circle", "Gauss", "rectangle", "exp_decay", "bartlett", "mfilter")[1],
  zoi_hl_ratio = 4,
  half_life = NULL,
  exp_decay_parms = c(1, 0.01),
  min_intensity = 0.01,
  max_dist = 50000,
  normalize = FALSE,
  extent_x_cut = terra::ext(x)[c(1,2)],
  extent_y_cut = terra::ext(x)[c(3,4)],
  na.policy = "omit",
  na.rm = TRUE,
  quiet = FALSE,
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

  # check if the input raster presents only a single value (1,NA)
  # if so, transform it into a binary map (1,0)
  r0 <- x
  if(use_terra) {
    if(diff(c(r0@ptr$range_min, r0@ptr$range_max)) == 0) r0 <- terra::classify(r0, cbind(NA, 0)) # binary map
  } else {
    if(diff(c(maxValue(r0), minValue(r0))) == 0) r0 <- raster::reclassify(r0, cbind(NA, 0)) # binary map
  }
  # plot(r0)

  # define filters
  if(type %in% c("exp_decay", "bartlett", "circle", "threshold", "step", "rectangle")) {
    if(length(zoi) == 1) {
      filt <- create_filter(r0, zoi = zoi, method = type,
                            zoi_hl_ratio = zoi_hl_ratio,
                            half_life = half_life,
                            max_dist = max_dist,
                            min_intensity = min_intensity,
                            normalize = normalize, ...)
    } else {
      filt <- purrr::map(zoi, function(z, ...) {
        create_filter(r0, zoi = z, method = type,
                      zoi_hl_ratio = zoi_hl_ratio,
                      half_life = half_life,
                      max_dist = max_dist,
                      min_intensity = min_intensity,
                      normalize = normalize, ...)
        })
    }
  }

  if(type == "mfilter") {
    filt <- zoi
  }

  # those methods were put above with the create_filter function
  # type = c("circle", "rectangle", "threshold", "step")
  # only Gauss is kept here, so far
  # if(type %in% c("circle", "Gauss", "rectangle", "threshold", "step")) {
  if(type %in% c("Gauss")) {

    # if(type %in% c("threshold", "step")) type <- "circle"
    # if(type == "rectangle") zoi = 2*zoi # for this case d is the side of the square

    if(length(zoi) == 1) {
      filt <- terra::focalMat(r0, d = zoi, type = type)
      if(!normalize) filt <- filt/max(filt, na.rm = T)
    } else {
      filt <- purrr::map(zoi, function(z) {
        ft <- terra::focalMat(r0, d = z, type = type)
        if(!normalize) ft <- ft/max(ft, na.rm = T)
        ft
      })
    }
  }

  # neighborhood analysis
  if(type == "mfilter") {
    # more than one matrix
    if("list" %in% class(filt)) {
      cuminf <- purrr::map2(filt, 1:length(zoi), function(f, z) {
        if(!quiet) print(paste0("Calculating for ZoI n. ", z, "..."))
        terra::focal(r0, w = f, na.policy = na.policy, na.rm = na.rm, ...)
      })
      if(use_terra) cumulative_r <- do.call(c, cuminf) else
        cumulative_r <- raster::stack(cuminf)
    } else {
      #only one matrix
      cumulative_r <- terra::focal(r0, w = filt, na.policy = na.policy, na.rm = na.rm, ...)
    }
  } else {
    if(length(zoi) == 1) {
      cumulative_r <- terra::focal(r0, w = filt, na.policy = na.policy, na.rm = na.rm, ...)
    } else {
      cuminf <- purrr::map2(filt, zoi, function(f, z) {
        if(!quiet) print(paste0("Calculating for ZoI = ", z, "..."))
        terra::focal(r0, w = f, na.policy = na.policy, na.rm = na.rm, ...)
      })
      if(use_terra) cumulative_r <- do.call(c, cuminf) else
        cumulative_r <- raster::stack(cuminf)
    }
  }

  # rename cumulative influence layer
  if(type == "mfilter") {
    if(!is.list(zoi)) name <- "influence_cumulative" else
      name <- paste0("influence_cumulative", 1:length(zoi))
  } else {
    name <- paste0("influence_cumulative_", type, zoi)
  }

  names(cumulative_r) <- name
  # should the result be plotted?
  if(plotit) plot(cumulative_r)

  # return cropped raster
  if(use_terra)
    terra::crop(cumulative_r, terra::ext(c(extent_x_cut, extent_y_cut)))
  else
    raster::crop(cumulative_r, raster::extent(c(extent_x_cut, extent_y_cut)))
}

calc_influence_cumulative_GRASS <- function(
  x,
  zoi = 100,
  type = c("circle", "Gauss", "rectangle", "exp_decay", "bartlett", "mfilter")[1],
  module = c("r.mfilter", "r.resamp.filter", "r.neighbors")[1],
  zoi_hl_ratio = 4,
  half_life = NULL,
  exp_decay_parms = c(1, 0.01),
  # hnorm_decay_parms = c(1, 20),
  min_intensity = 0.01,
  max_dist = 50000,
  divisor = 1,
  normalize = FALSE,
  extent_x_cut = NULL,
  extent_y_cut = NULL,
  parallel = TRUE,
  output_map_name = NULL,
  remove_intermediate = TRUE,
  overwrite = FALSE,
  quiet = TRUE,
  ...) {

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

  ##### CUT extent to be implemented

  # start by setting the region
  rgrass7::execGRASS("g.region", raster = x, flags = flags_region)

  # check if the input raster presents only a single value (1,NA)
  # if so, transform it into a binary map (1,0) or (integer numbers, 0)
  input_bin <- x
  values_input <- rgrass7::execGRASS("r.category", map = input_bin, separator = " ", flags = "quiet")
  values_input <- sort(as.numeric(attributes(values_input)$resOut))
  if(sort(values_input) != c(0,1)) {
    stop("Please make sure the input map is binary, with values 0 or 1 only.")
    # input_bin <- paste0(x, "_bin")
    # if(overwrite_bin) {
    #   fl_bin <- ifelse("overwrite" %in% flags, flags, c(flags, "overwrite"))
    #   rgrass7::execGRASS("r.mapcalc", expression = sprintf("%s = %s", input_bin, x), flags = fl_bin)
    # }
    # rgrass7::execGRASS("r.null", map = input_bin, null = 0)
    # if(remove_intermediate) to_remove <- c(to_remove, input_bin)
  }

  # define the name of the output map
  if(!is.null(output_map_name)) {
    # given name if this is given as a parameter
    out_map <- output_map_name
  } else {
    # define name as input + cumulative + method
    out_map = paste0(x, "_inf_cumulative_", type)

    # add zoi
    # if(type != "mfilter") {
    #   # single or multiple strings, if zoi is 1 or nultiple values
    #   out_map <- paste0(out_map, zoi)
    # }
  }

  # get resolution
  region <- rgrass7::gmeta()
  resolution <- region$nsres

  # perform calculations for "r.resamp.filter"
  # if(module == "r.resamp.filter") {
  #   if(type == )
  # }

  # for r.neighbors, it always performs the average, it is not possible to sum
  # the matrix is always normalized
  # one must use the argument size in conjunction with the weight matrix
  # we must a specific output from save_mfilter for that, with only the numbers
  # r.neighbors input=private_cabins_sub_bin output=test_neighbors method=count size=21 weight=test_neighbors_exp_filt500.txt --o

  # perform calculations for "r.mfilter"
  if(module == "r.mfilter") {

    # define filters
    if(type %in% c("exp_decay", "bartlett", "circle", "threshold", "step", "rectangle")) {
      filter_count <- zoi
      filter_file <- tempfile(paste0("my_filter", filter_count, "_"))
      if(length(zoi) == 1) {
        filt <- create_filter(r = resolution, zoi = zoi, method = type,
                              zoi_hl_ratio = zoi_hl_ratio,
                              half_life = half_life,
                              max_dist = max_dist,
                              min_intensity = min_intensity,
                              divisor = divisor,
                              normalize = normalize, save_txt = TRUE,
                              save_file = filter_file, ...)
      } else {
        filt <- purrr::map2(zoi, filter_file, function(z, file, ...) {
          create_filter(r = resolution, zoi = z, method = type,
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

    # In these cases the matrix is not defined by create_filter
    if(type %in% c("mfilter", "Gauss")) {

      # Filters pre-defined for "mfilter"
      if(type == "mfilter") {
        # set
        filt <- zoi
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

      # create filters with focalMat for "Gauss", zoi represents the sd
      if(type == "Gauss") {

        # create raster to define resolution
        reg_proj <- rgrass7::getLocationProj()
        r0 <- terra::rast(nrows = region$rows, ncols = region$cols,
                          xmin = region$w, xmax = region$e,
                          ymin = region$s, ymax = region$n,
                          crs = reg_proj)

        if(length(zoi) == 1) {
          filt <- terra::focalMat(r0, d = zoi, type = type)
          if(!normalize) filt <- filt/max(filt, na.rm = T)
        } else {
          filt <- purrr::map(zoi, function(z) {
            ft <- terra::focalMat(r0, d = z, type = type)
            if(!normalize) ft <- ft/max(ft, na.rm = T)
            ft
          })
        }
      }

      # matrices created, save them outside R for use in GRASS GIS

      # for one matrix only
      if(is.matrix(filt)) {
        filter_count <- 1
        filter_file <- tempfile(paste0("my_filter_", type, filter_count, "_"))
        # save matrix outside R for use within GRASS GIS
        save_mfilter(filt, zoi = "", method = type,
                     divisor = divisor,
                     save_format = c("GRASS_r.mfilter"),
                     save_file = filter_file,
                     parallel = parallel,
                     separator = " ")
      } else {
        # for multiple matrices
        if(is.list(filt)) {
          filter_count <- 1:length(filt)
          filter_file <- tempfile(paste0("my_filter_", type, filter_count, "_"))
          # save matrices outside R for use within GRASS GIS
          filt <- purrr::map2(filt, filter_file, function(f, file, ...) {
            save_mfilter(f, zoi = "", method = type,
                         divisor = divisor,
                         save_format = c("GRASS_r.mfilter"),
                         save_file = file,
                         parallel = parallel,
                         separator = " ")
          })
        }
      }
    }

    # set parameters for neighborhood analysis

    # only one matrix
    if(is.matrix(filt)) {
      out_names <- out_map
      if(type != "mfilter") out_names <- paste0(out_names, zoi)
      parms <- list(input = input_bin, output = out_names, filter = filter_file)
    } else {
      # several matrices or zoi values
      if(is.list(filt)) {
        parms <- purrr::map2(filter_count, filter_file, function(x, y)
          list(input = input_bin, output = paste0(out_map, x), filter = y))
        out_names <- purrr::map(parms, ~ .$output) %>% unlist()
      }
    }
  }

  # run neighborhood analysis
  if("list" %in% class(filt)) {

    # loop for matrices or zoi values
    for(i in 1:length(filt)) {
      parm <- parms[[i]]
      z <- filter_count[[i]]

      # message
      if(!quiet) print(paste0("Calculating for ZoI ", z, "..."))
      # region
      # set region
      rgrass7::execGRASS("g.region", raster = parm$input, flags = flags_region)
      # calculate
      rgrass7::execGRASS(module, parameters = parm, flags = flags)
    }

  } else {

    # for only one matrix/zoi value
    parm <- parms
    z <- filter_count

    # message
    if(!quiet) print(paste0("Calculating for ZoI ", z, "..."))
    # region
    # set region
    rgrass7::execGRASS("g.region", raster = parm$input, flags = flags_region)
    # calculate
    rgrass7::execGRASS(module, parameters = parm, flags = flags)

  }

  # return only names
  return(out_names)
}
