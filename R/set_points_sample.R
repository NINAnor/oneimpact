#' Simulate regular or random points in 2D
#'
#' This function simulates the coordinates of points either regularly or
#' randomly  distributed in two dimensions. Other point patterns 
#' can be also chosen, see the `type` argument for the [sf::st_sample]
#' function. The points are generated within an input polygon or the
#' bounding box of an input raster. If no spatial object is used as input,
#' The size of the landscape is defined by the `extent_x` and `extent_y` parameters. 
#' It assumes the landscape is square or rectangularly shaped.
#' 
#' @param n_features `[integer(1)=1000]` \cr Total number of points to spread in space.
#' @param type `[character(1)="regular"]{"regular", "random"}` Pattern for the creation
#' of points is space. Other methods are also accepted, check the `type` argument 
#' for the [sf::st_sample] function.
#' @param base_polygon `[RasterLayer,sfc_POLYGON]` \cr Polygon (from class `sf` or `sfc`)
#' inside which the points will be created. If a `RasterLayer`, the `bbox` of the 
#' raster is used as a polygon.
#' @param extent_x,entent_y `[numeric vector(2)=c(0,1)]` \cr Vector representing the minimum and
#' maximum extent in x and y within which the points should be placed, in the format c(min,max).
#' 
#' @return A `data.frame` with the (x,y) random coordinates.
#' 
#' @examples 
#' library(sf)
#' 
#' pts <- set_points_sample(100)
#' plot(pts)
#' pts2 <- set_points_sample(100, type = "random")
#' plot(pts2)
#' 
#' library(terra)
#' library(stars)
#' x <- rast(system.file("external/test.grd", package="raster"))
#' pts3 <- set_points_sample(100, base_polygon = x)
#' plot(pts3)
#' 
#' @export
set_points_sample <- function(n_features = 1000,
                              type = c("regular", "random")[1],
                              base_polygon = NULL,
                              extent_x = c(0,1),
                              extent_y = extent_x) {
  
  if(is.null(base_polygon)) {
    # Assumes a square/rectangular shape and creates a polygon 
    coord <- expand.grid(x = extent_x, y = extent_y)
    coord <- rbind(coord[1:2,], coord[4,], coord[3,], coord[1,]) %>% as.matrix()
    pol <- sf::st_polygon(list(coord))
    
  } else {
    # Creates bbox around a raster
    if(class(base_polygon) %in% c(paste0("Raster", c("Layer", "Stack", "Brick")), "SpatRaster") ) {
      pol <- base_polygon %>% 
        sf::st_bbox() %>% 
        sf::st_as_sfc() 
    } else {
      pol <- base_polygon
    }
    
  }
  
  # sample points
  pts <- sf::st_sample(pol, size = n_features, type = type)
  
  # return coordinates
  data.frame(st_coordinates(pts))
}
