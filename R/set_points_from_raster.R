#' Simulate points using input raster as weights 
#'
#' This function simulates point patterns in space using the values of
#' an input raster as weights or probabilities for selecting a point in a 
#' given location. It was designed to simulate points based on  neutral landscape 
#' models but it works with other input rasters as well.
#' 
#' The function works by first selecting random pixels in the landscape and 
#' finding their centers, then adding random variation within each pixel to
#' define the final point locations. 
#' It was based on this StackExchange very useful answer from "Spacedman":
#' https://gis.stackexchange.com/questions/224321/randomly-generate-points-using-weights-from-raster
#'
#' TO IMPROVE: implement with terra package
#'
#' @param n_features `[integer(1)=1000]` \cr Total number of features to spread in space.
#' @param base_raster `[RasterLayer]` \cr Input raster used for defining the weights.
#'
#' @returns The coordinates (x,y) of the simulated points.
#'
#' @example examples/set_points_from_raster_example.R
#'
#' @export

# function to simulate points using input raster as weights 
set_points_from_raster <- function(base_raster, n_features = 1000) {
  
  # get parameters
  res <- terra::res(base_raster)
  
  # random points in the center of the cells
  ptscell <- sample(1:ncell(base_raster), n_features, prob = base_raster[], replace = TRUE)
  # get the centers
  center <- xyFromCell(base_raster, ptscell)
  # add random values within the pixels
  pts <- center + cbind(runif(nrow(center), - res[1]/2, res[1]/2),
                        runif(nrow(center), - res[2]/2, res[2]/2))
  
  # return the points
  data.frame(pts)
}
