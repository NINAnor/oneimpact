#' Preparing data for spatially stratified cross‐validation schemes
#'
#' This function sets up the variables used spatial stratification sampling to
#' be used in the bootstrapped modeling approach in [oneimpact].
#' The function defines the hierarchical levels H0,
#' used for block cross-validation, as well as the levels H1 and H2, used for setting the
#' spatial structure of the sampling for model fitting and tuning.
#' Levels H1 and and H2 are organized in hierarchical spatial blocks. Blocks H1 correspond
#' to larger square blocks, each of which is subdivided into four, smaller H2 blocks
#' (with side equal `block_size`). H2 blocks are separated
#' spatially so that 3 of them are used for model fitting the the fourth is used for
#' model tuning.
#' The level H0 comes from the data set and should represent some other level for
#' block cross validation, such as population ID, animal ID, year or study area.
#' The block H0 is used for drawing validation sets, to allow a thorough evaluation
#' of the fitted models in all blocks H0.
#'
#' The function returns a `data.frame` with the blocks for the hierarchical
#' levels H0, H1, and H2 for each observation in the input data set,
#' to be used to create samples for the bootstrapped modeling
#' approach using the [oneimpact::create_samples()] function, in case spatial
#' stratified samples are desired.
#'
#' @note
#' To be implemented for input = data.frame
#' #terra::vect(data[data$case==1,], geom = c("x", "y"))
#' To be implemented for track objects - already have crs
#' Put H0 here as well.
#'
#' @param x `[data.frame,sf,SpatVector]` \cr Vector of points to be spatially stratified, as a
#' [sf] or a [terra::SpatVector] object. If a `data.frame`, the columns corresponding to the (x,y) coordinates must be given in
#' `coords`.
#' @param colH0 `[numeric,character=NULL]` \cr Column number or name to define the ids of the H0 level
#' - the one with ecological meaning, e.g. individuals, populations, or study areas, used for block cross-validating
#' the predictions of the fitted models. Default is `NULL`, in which case there is no block H0 defined.
#' @param col_id `[numeric,character=NULL]` Column number or name with the ID of the rows of the
#' data observations. In step-selection analysis, this should be the column showing the
#' number of the strata of each step. For resource selection analysis and environmental niche modeling,
#' this might be the row id, for instance.
#' @param H1_as_H0 `[logical(1)=FALSE]` \cr Whether the spatial blocks of level H1 should be used
#' as the block H0, in case no block H0 is provided (if `colH0 = NULL`). This parameter is ignored
#' if `colH0` is provided.
#' @param k `[numeric(1)=4]` \cr Number of H2 blocks within each block H1. Should be 4, 9, 16, or some
#' number `k = x**2`, where x is an integer > 1. Default is `k = 4`.
#' TO BE IMPLEMENTED: Number of parts for k-fold cross validation within H1 hierarchical level, for
#' tuning (setting the penalty parameter). Could be used for nested cross-validation.
#' @param block_size `[numeric(1)=10000]` \cr Size (side of a square) of the blocks for H2 level,
#' map units (generally meters). The size of the H1 level blocks is defined as `sqrt(k)*block_size`.
#' Default is `block_size = 10000`
#' @param buffer `[numeric(1)=1000]` \cr Buffer added around the points before creating the blocks,
#' to make sure all points are included in the samples. Default is `buffer = 1000`.
#' @param coords [string,vector] \cr Vector with the names of the columns with the (x,y)
#' coordinates of the locations from the data set. Default is `NULL`, in which case `x`
#' should be a `sf` or a `terra::SpatVector` object. If `x` is a `data.frame`, `coords`
#' must be provided.
#' @param all_cols `[logical=FALSE]` \cr If `TRUE`, and if `x` is a `data.frame`,
#' the spatial strata blocks are appended as columns in the input data `x`.
#' @param crs `[string=""]` \cr Coordinate reference system (CRS) of the observations,
#' if `x` is a `data.frame`.
#' @param plot_grid `[logical=TRUE]` \cr if `TRUE` (default), the grid with spatial blocks and
#' observations is plotted.
#' @param save_grid `[character=NA]{NA, "raster", "vector"}` \cr Should the grid which defines
#' the H1 and H2 blocks be saved? NOT IMPLEMENTED.
#'
#' @return A `data.frame` with the blocks at hierarchical levels H0, H1, and H2 corresponding to
#' each of the observations in the input data set `x`.
#'
#' @examples
#' data(reindeer)
#' library(terra)
#' library(amt)
#'
#' # no block H0, spatial blocks
#' spat_strat(reindeer, block_size = 5000, coords = c("x", "y"))
#'
#' # with block H0
#' spst <- spat_strat(reindeer, coords = c("x", "y"), colH0 = "original_animal_id",
#'                    all_cols = TRUE)
#' # Visualize level H0 - individuals
#' spst_vect <- terra::vect(spst, geom = c("x", "y"))
#' terra::plot(spst_vect, "blockH0")
#' # Visualize level H1
#' terra::plot(spst_vect, "blockH1", type = "classes")
#' # Visualize level H2 for blockH1 numbers 6 to 10
#' terra::plot(spst_vect, col = grey(0.7))
#' terra::plot(spst_vect[spst$blockH1 == 6], "blockH2", type = "classes", add = TRUE) # only 6
#' terra::plot(spst_vect[spst$blockH1 %in% c(6,7,10,11)], "blockH2", type = "classes", add = TRUE) # 6-10
#' terra::plot(spst_vect[spst$blockH1 == 6], "blockH2", type = "classes") # zoom to 6
#'
#' @export
spat_strat <- function(x,
                       colH0 = NULL,
                       colID = NULL,
                       H1_as_H0 = FALSE,
                       k = 4,
                       block_size = 10000,
                       buffer = 1000,
                       coords = NULL,
                       all_cols = FALSE,
                       crs = "",
                       plot_grid = TRUE,
                       save_grid = c(NA_character_, "raster", "vector")[1]) {

  # check input class for multiple uses
  if(all(class(x) != "SpatVector")) {
    if(any(class(x) %in% c("sf", "SpatialPoints", "SpatialPointsDataFrame"))) {
      # convert x to SpatVector
      y <- terra::vect(x)
    } else {
      if(any(class(x) %in% c("data.frame"))) {
        if(is.null(coords)) {
          stop("For `data.frame` objects, you must specify the (x,y) columns in the `coords` argument.")
        } else {
          y <- terra::vect(as.data.frame(x), geom = coords, crs = crs)
        }
      }
    }
  }

  # set points to blocks of levels H1 and H2
  blocks <- spat_strat_blocks(y, k = k, block_size = block_size, max_sl = buffer,
                              plot_grid = plot_grid, save_grid = save_grid)

  # set blockH0 according to column colH0 and class of input x
  if(is.null(colH0)) {

    if(H1_as_H0) {
      blockH0 <- blocks[["blockH1"]]
    } else {
      blockH0 <- NA_integer_
    }

  } else{

    if(all(class(x) != "SpatVector")) {
      blockH0  <- x[[colH0]]
    } else {
      blockH0  <- x[[colH0]][,1]
    }
  }

  # put block H0, H1, and H2 together
  blocks <- cbind(blocks, data.frame(blockH0 = blockH0))
  sp_strat <- blocks[, c("id", "blockH0", "blockH1", "blockH2")]

  if(!is.null(colID)) {
    if(all(class(x) != "SpatVector")) {
      sp_strat$id  <- x[[colID]]
    } else {
      sp_strat$id  <- x[[colID]][,1]
    }
  }

  # add to original data.frame
  if(any(class(x) %in% c("data.frame")) & all_cols) {
    sp_strat <- cbind(x, sp_strat)
  }

  sp_strat
}

spat_strat_blocks <- function(x,
                              k = 4,
                              block_size = 10000,
                              max_sl = 1000,
                              plot_grid = T,
                              save_grid = c(NA_character_, "raster", "vector")[1]) {

  # get extent of x
  ext_orig <- terra::ext(x)

  # set k
  k_use <- round(sqrt(k))**2 # round since we want quadrants
  if(k_use == 1) k_use <- k_use <- 4
  if(k_use != k) warning(paste0("Impossible to use k = ", k, " in a spatial setup. Using k = ", k_use, " instead."))
  aggr_factor <- sqrt(k_use)

  # round block size
  block_size <- 2 * round(block_size/2)

  # define xmin, ymin, and range, rounding to make sure to encompass all
  # start and end points and steps
  xlimit <- unname(range(ext_orig)[1] + 2*max_sl)
  ylimit <- unname(range(ext_orig)[2] + 2*max_sl)
  xlimit <- ((xlimit %/% block_size)+1)*block_size
  ylimit <- ((ylimit %/% block_size)+1)*block_size
  rounding <- 10^(floor(log10(block_size))-1)
  xmin <- floor((min(ext_orig)[1] - max_sl)/rounding)*rounding
  ymin <- floor((min(ext_orig)[2] - max_sl)/rounding)*rounding

  # Define hierarchical level H2
  ext_out <- terra::ext(c(xmin, xmin+xlimit, ymin, ymin+ylimit))
  xcols <- xlimit/(block_size/2)
  yrows <- ylimit/(block_size/2)

  H2rast <- terra::rast(ext_out, nrows = yrows, ncols = xcols, crs = terra::crs(x))
  terra::values(H2rast) <- 1:(nrow(H2rast)*ncol(H2rast))
  #plot(H2rast)
  #plot(x, add=T)

  H1rast <- terra::aggregate(H2rast, fact = aggr_factor, fun = mean)
  terra::values(H1rast) <- 1:(nrow(H1rast)*ncol(H1rast))
  #plot(H2rast)
  #plot(H1rast)

  ptsH1 <- terra::extract(H1rast, x)
  ptsH2 <- terra::extract(H2rast, x)
  blocks <- cbind(ptsH1, ptsH2[,2])
  names(blocks) <- c("id", "blockH1", "blockH2")

  if (plot_grid){
    terra::plot(H1rast)
    terra::plot(x, add = T)
  }

  return(blocks)
}

# Spatially stratified cross‐validation is less likely than randomly partitioned cross‐validation to
# overestimate predictive performance
# due to spatial autocorrelation (Radosavljevic & Anderson, 2014; Veloz, 2009; Wenger & Olden, 2012),
# so it is a good way to assess predictive power when independent evaluation data are not available.


# test subsetting different objects
# data(reindeer)
# col <- "sex"
# # data.frame
# reindeer[[col]]
# # sf
# library(sf)
# sf::st_as_sf(reindeer, coords = c("x", "y"), remove = FALSE)[[col]]
# # terra
# library(terra)
# terra::vect(reindeer, geom = c("x", "y"))[[col]][,1]
# # sp
# library(sp)
# SpatialPointsDataFrame(reindeer[,c("x", "y")], reindeer)[[col]]
