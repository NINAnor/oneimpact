# Isolation and mean isolation of points in space

Measures of isolation and mean isolation to a set of points in space.
`isolation()` creates random points in a landscape and calculates the
nearest neighbor distance from each of them to another set of points
passed as input, `x`. `mean_isolation()` calculates the average
isolation calculated through `isolation()`.

## Usage

``` r
isolation(x, n_rand = 100, ext = c(0, 1, 0, 1), lonlat = FALSE)

mean_isolation(x, n_rand = 100, ext = c(0, 1, 0, 1), lonlat = FALSE)
```

## Arguments

- x:

  `[data.frame]`  
  `data.frame` with (x,y) coordinates in the columns.

- n_rand:

  `[numeric(1)=100]`  
  Number of random points to be created in space, to compute the
  distance to `x`.

- ext:

  `[numeric(4)=c(0,1,0,1)]`  
  Extent of the space within which the random positions should be
  created: `c(xmin, xmax, ymin, ymax)`.

- lonlat:

  `[logical(1)=FALSE]`  
  Whether the distance between points should be calculated in an WGS
  ellipsoid (`lonlat = TRUE`) or on a plane (`lonlat = FALSE`). See
  [`raster::pointDistance()`](https://rdrr.io/pkg/raster/man/pointDistance.html)
  for more details.

## Value

`isolation()` returns the distance from each random point to the nearest
neighbor point in `x`. `mean_isolation()` returns the average nearest
neighbor distance from all random positions to the points in `x`.

## Details

So far the function only works for a square landscape. In the future we
can implement that for polygons or rasters with masks or null cells if
necessary, in an approach similar to set_points_sample.

## Examples

``` r
pts <- set_points(n_features = 100, method = "random", centers = 1, width = 0.1)[[1]]
isolation(pts)
#>   [1] 0.018089446 0.090221469 0.035366194 0.020597131 0.046765397 0.093244627
#>   [7] 0.106545038 0.105277314 0.055122624 0.030620770 0.046711020 0.042383009
#>  [13] 0.014562017 0.045221034 0.028529082 0.067636651 0.035080641 0.032090491
#>  [19] 0.041326978 0.100999713 0.042862374 0.026648799 0.070148659 0.064095964
#>  [25] 0.020064086 0.016014413 0.019934726 0.045951283 0.057659148 0.049261680
#>  [31] 0.062383647 0.026882708 0.037342477 0.032587149 0.032242097 0.069560223
#>  [37] 0.057012786 0.039282515 0.044351664 0.042618699 0.044839194 0.017875209
#>  [43] 0.085758998 0.051898315 0.042626919 0.050329569 0.015699021 0.034438381
#>  [49] 0.032308239 0.016297550 0.015776415 0.066926295 0.070365968 0.111154575
#>  [55] 0.098805194 0.036904638 0.061876999 0.068296818 0.097268087 0.065078225
#>  [61] 0.072724306 0.099610316 0.053759339 0.009697322 0.062353489 0.027793917
#>  [67] 0.014637550 0.029367169 0.051903128 0.011217081 0.072929554 0.061199390
#>  [73] 0.037235810 0.018697607 0.038175864 0.062960375 0.059988320 0.043891634
#>  [79] 0.015193774 0.020632509 0.041528786 0.049099589 0.061626925 0.089678810
#>  [85] 0.024178623 0.058766316 0.051812620 0.084649489 0.079431968 0.020324333
#>  [91] 0.071308027 0.058023803 0.015883813 0.048977570 0.032135436 0.105796959
#>  [97] 0.070599789 0.062221706 0.045613226 0.035624219 0.081454856
mean_isolation(pts)
#> [1] 0.05125018
```
