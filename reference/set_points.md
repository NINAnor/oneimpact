# Simulate points in a landscape

This function simulates point patterns in space and rasterize them. The
idea is to mimic the spatial distribution of point-type infrastructure,
such as houses, cabins, or turbines, for instance. The function returns
a list with the position of the points and a binary raster with 1 where
there are points and NA elsewhere. If created with a raster to define
the weights, this base raster is also returned in the output.

## Usage

``` r
set_points(
  n_features = 1000,
  method = c("mobsim", "regular", "random", "raster")[1],
  centers = 1,
  width = 0.05,
  base_raster = NULL,
  point_coordinates = NULL,
  res = 0.1,
  extent_x = c(0, 1),
  extent_y = c(0, 1),
  buffer_around = 0,
  return_base_raster = TRUE,
  use_terra = TRUE,
  crs = ""
)
```

## Arguments

- n_features:

  `[integer(1)=1000]`  
  Total number of features to spread in space.

- method:

  `[character(1)]{"mobsim", "regular", "random", "raster"}`  
  Method used to simulate points in space. `mobsim` uses the function
  [`mobsim::sim_thomas_community()`](https://rdrr.io/pkg/mobsim/man/sim_thomas_community.html)
  from the [mobsim](https://rdrr.io/pkg/mobsim/man/mobsim-package.html)
  package to simulate points. `raster` uses a base raster map as input
  to define weights and simulate the random points.

- centers:

  `[integer(1)=1]`  
  Number of centers around which the features will be placed. Used only
  if `method = "mobsim"`.

- width:

  `[numeric(1)=0.05]`  
  Mean distance between each of the features in a cluster and the center
  of the cluster. Used only if `method = "mobsim"`.

- base_raster:

  `[RasterLayer=NULL]`  
  Base raster to define weights for creating the random points. Used
  only if `method = "raster"`.

- point_coordinates:

  `[data.frame=NULL]`  
  `data.frame` with (x,y) columns with coordinates already taken from
  elsewhere. This option is intended for when the points' coordinates
  were already generated or taken from a real landscape. In this case,
  no points are simulated and they are just rasterized (so that
  distances or other derived variables might be calculated).

- res:

  `[numeric(1)=0.1]`  
  Resolution of the output raster.

- extent_x, entent_y:

  `[numeric vector(2)=c(0,1)]`  
  Vector representing the minimum and maximum extent in x and y within
  which the points should be placed, in the format c(min,max).

- buffer_around:

  `[numeric(1)=0.1]`  
  Size of the buffer around the extent of the landscape, to avoid edge
  effects when calculating densities using neighborhood analysis.

- return_base_raster:

  `[logical(1)=TRUE]`  
  Whether the base_raster should be returned in the output list. This is
  `NULL` for `method = "mobsim"`.

- use_terra:

  `[logical(1)=TRUE]`  
  If `TRUE` (default), the `rast` element created from the points is a
  `SpatRaster` object from `terra` package is created. If `FALSE`, it is
  a `RasterLayer` from `raster` package is created.

- crs:

  `[character(1)]`  
  Specification for the coordinate reference system of the `rast` object
  created from the points. Default is
  `"+proj=utm +zone=1 +datum=WGS84"`. argument.

## Value

A list with three elements: (1) `pts`, the coordinates (x,y) of the
simulated points; (2) `rast`, a binary raster containing the landscape,
with 1 where there points and NA elsewhere; (3) `base_rast`, the base
raster used to weigh the simulation of points. If `method = "mobsim"` or
`"regular"` or `"random"`, `base_rast` is `NULL`.

## Details

If `method = "mobsim"`, the function builds upon the function
[`mobsim::sim_thomas_community()`](https://rdrr.io/pkg/mobsim/man/sim_thomas_community.html)
from the [mobsim](https://rdrr.io/pkg/mobsim/man/mobsim-package.html)
package. Originally the function is intended to simulate positions of
multiple species in the context of species abundance distribution
studies, but it fits well in case of a single species (or point patterns
for a single type of feature). In this case, the points are simulated
based on the number of centers/patches of points and their width.

If `method = "raster"`, the function uses an input raster (defined by
the argument `base_raster`) to define the probabilities of setting a
given point in a certain pixel in space.

## Examples

``` r
#-----
# using mobsim
library(terra)
library(mobsim)
#> 
#> Attaching package: ‘mobsim’
#> The following object is masked from ‘package:oneimpact’:
#> 
#>     dist_decay

set.seed(1234)

# gradient distribution
ext <- 30000
wd <- ext/5
pts <- set_points(n_features = 1000, centers = 1,
                  width = wd, res = 100,
                  extent_x = c(0, ext), extent_y = c(0, ext))
plot(pts$pts)

plot(pts$rast, col = "black")


# one focus of features, with buffer around
wd <- ext/20
pts <- set_points(n_features = 1000, centers = 1,
                  width = wd, res = 100,
                  extent_x = c(0, ext), extent_y = c(0, ext),
                  buffer_around = 10000)
plot(pts$pts)

plot(pts$rast, col = "black")


#-----
# using base raster

# raster
set.seed(12)
r <- raster::raster(matrix(runif(12),3,4)) |>
  raster::disaggregate(fact = 10)

# points from raster
pts <- set_points(n_features = 100, method = "raster",
                  base_raster = r)
plot(pts$base_rast)

plot(pts$pts)

plot(pts$rast, col = "black")


#-----
# using random or regular

set.seed(123)
ext <- 30000
pts <- set_points(n_features = 1000, method = "random",
                  res = 100,
                  extent_x = c(0, ext), extent_y = c(0, ext))
plot(pts$pts)

plot(pts$rast, col = "black")


pts <- set_points(n_features = 1000, method = "regular",
                  res = 100,
                  extent_x = c(0, ext), extent_y = c(0, ext))
plot(pts$pts)

plot(pts$rast, col = "black")


#-----
# using point coordinates as input
pt_input <- data.frame(x = c(0.5, 0.7), y = c(0.5, 0.3))
pts <- set_points(point_coordinates = pt_input)
plot(pts$pts)

plot(pts$rast, col = "black")
```
