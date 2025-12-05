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

  `[numeric(x)=c(0,1)]` Extent of the space within which the random
  positions should be created c(x or ymin, x or ymax).

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
#>   [1] 0.082082527 0.086352651 0.092660941 0.033989324 0.050086774 0.044043924
#>   [7] 0.079940197 0.069154057 0.012799689 0.088310023 0.133969121 0.052798728
#>  [13] 0.064125990 0.061266780 0.033580613 0.024674247 0.052281233 0.066475139
#>  [19] 0.034947389 0.080351889 0.077491703 0.065512901 0.070029064 0.035174023
#>  [25] 0.053940119 0.100590540 0.072440480 0.071537488 0.079876374 0.090779646
#>  [31] 0.039714024 0.054202575 0.020691884 0.044572917 0.026994918 0.027156772
#>  [37] 0.031950903 0.025883465 0.030080789 0.070052510 0.044753008 0.049006889
#>  [43] 0.049203603 0.018808723 0.028172559 0.049025817 0.048342019 0.080413225
#>  [49] 0.055897993 0.053068308 0.035125866 0.021457389 0.057871282 0.077686753
#>  [55] 0.015040050 0.039613286 0.039447330 0.030826041 0.056072132 0.024553416
#>  [61] 0.035314345 0.090709131 0.062038550 0.043226221 0.056666972 0.030510302
#>  [67] 0.007637838 0.111346564 0.072305701 0.056037034 0.038055308 0.062466313
#>  [73] 0.013840698 0.012342582 0.050456630 0.051113391 0.066767501 0.019358712
#>  [79] 0.042892294 0.030431971 0.071089956 0.088874554 0.064598630 0.084962641
#>  [85] 0.035952402 0.013454626 0.096874271 0.037617063 0.071481711 0.025236252
#>  [91] 0.058576248 0.070009815 0.062291998 0.066426509 0.071089958 0.028168172
#>  [97] 0.063957029 0.034056139 0.021865841 0.053348144 0.097344844
mean_isolation(pts)
#> [1] 0.05343212
```
