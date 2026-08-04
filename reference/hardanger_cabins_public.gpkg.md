# Public cabins vector data for the Hardangervidda wild reindeer area in Norway

Dataset containing the location of large, public DNT cabins in the
surroundings of the Hardangervidda wild reindeer area. Retrieved from
the public N50 dataset.

## Format

A geopackage file. Projected CRS: ETRS89 / UTM zone 33N. The vector file
presents the following columns:

- gid: Line number, corresponding to the original dataset

- buildtype: Type of building (code) in the original dataset

- city: Code of the municipality where the cabin is located

- value: Value 1, to be used for rasterization purposes

## Source

<https://register.geonorge.no/det-offentlige-kartgrunnlaget/n50-kartdata/ea192681-d039-42ec-b1bc-f3ce04c189ac>

## Examples

``` r
(f <- system.file("vector/hardanger_cabins_public.gpkg", package = "oneimpact"))
#> [1] "/home/runner/work/_temp/Library/oneimpact/vector/hardanger_cabins_public.gpkg"
v <- sf::st_read(f)
#> Reading layer `hardanger_cabins_public' from data source 
#>   `/home/runner/work/_temp/Library/oneimpact/vector/hardanger_cabins_public.gpkg' 
#>   using driver `GPKG'
#> Simple feature collection with 21 features and 4 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 63357.86 ymin: 6675293 xmax: 132350.9 ymax: 6726508
#> Projected CRS: ETRS89 / UTM zone 33N
plot(v[1])

```
