# Simulate regular or random points in 2D

This function simulates the coordinates of points either regularly or
randomly distributed in two dimensions. Other point patterns can be also
chosen, see the `type` argument for the
[`sf::st_sample()`](https://r-spatial.github.io/sf/reference/st_sample.html)
function. The points are generated within an input polygon or the
bounding box of an input raster. If no spatial object is used as input,
The size of the landscape is defined by the `extent_x` and `extent_y`
parameters. It assumes the landscape is square or rectangularly shaped.

## Usage

``` r
set_points_sample(
  n_features = 1000,
  type = c("regular", "random")[1],
  base_polygon = NULL,
  extent_x = c(0, 1),
  extent_y = extent_x
)
```

## Arguments

- n_features:

  `[integer(1)=1000]`  
  Total number of points to spread in space.

- type:

  `[character(1)="regular"]{"regular", "random"}` Pattern for the
  creation of points is space. Other methods are also accepted, check
  the `type` argument for the
  [`sf::st_sample()`](https://r-spatial.github.io/sf/reference/st_sample.html)
  function.

- base_polygon:

  `[RasterLayer,sfc_POLYGON]`  
  Polygon (from class `sf` or `sfc`) inside which the points will be
  created. If a `RasterLayer`, the `bbox` of the raster is used as a
  polygon.

- extent_x, entent_y:

  `[numeric vector(2)=c(0,1)]`  
  Vector representing the minimum and maximum extent in x and y within
  which the points should be placed, in the format c(min,max).

## Value

A `data.frame` with the (x,y) coordinates of the simulated points.

## Examples

``` r
library(sf)

pts <- set_points_sample(100)
plot(pts)

pts2 <- set_points_sample(100, type = "random")
plot(pts2)


library(raster)
x <- rast(system.file("external/test.grd", package="raster"))
pts3 <- set_points_sample(100, base_polygon = x)
plot(pts3)

```
