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
#>   [1] 0.043699719 0.062651997 0.083621078 0.029002141 0.056218695 0.058303430
#>   [7] 0.103153073 0.006549776 0.038288423 0.028468476 0.083004952 0.103426887
#>  [13] 0.015947344 0.015851437 0.129169666 0.017296323 0.086540345 0.091112309
#>  [19] 0.021166060 0.037133030 0.099901574 0.075795043 0.099501910 0.110305055
#>  [25] 0.036671237 0.056665739 0.035108415 0.035834801 0.057468041 0.045191709
#>  [31] 0.037145431 0.038703942 0.036184402 0.037399072 0.044781946 0.051535996
#>  [37] 0.023328058 0.088453606 0.078246381 0.085317892 0.024047993 0.056814385
#>  [43] 0.063855475 0.065902219 0.055728898 0.040573150 0.039561394 0.062923265
#>  [49] 0.053880631 0.074799998 0.075406580 0.107228992 0.070988682 0.107201657
#>  [55] 0.017460893 0.115730825 0.035595276 0.057171342 0.052013942 0.009176495
#>  [61] 0.065285040 0.072698105 0.036982613 0.042598787 0.091264864 0.033039342
#>  [67] 0.077017389 0.049300816 0.075652621 0.082807113 0.066946510 0.090167922
#>  [73] 0.064818280 0.030299529 0.026436968 0.128588775 0.030387029 0.073152063
#>  [79] 0.002915024 0.018614491 0.055157651 0.046825627 0.049347757 0.100850900
#>  [85] 0.161595151 0.038126475 0.083522809 0.060266137 0.049997637 0.058276982
#>  [91] 0.067871533 0.042676467 0.030062931 0.008632426 0.036947197 0.102735835
#>  [97] 0.093800601 0.027836436 0.013944608 0.071492978 0.109577525
mean_isolation(pts)
#> [1] 0.04880009
```
