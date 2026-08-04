# Tourist trails vector data for the Hardangervidda wild reindeer area in Norway

Dataset containing the location of tourist trails in the surroundings of
the Hardangervidda wild reindeer area. Retrieved from the public N50
dataset.

## Format

A geopackage file. Projected CRS: ETRS89 / UTM zone 33N. The vector file
presents the following columns:

- gid: Line number/identifier

- area: Name of the wild reindeer area, if within one

- traffic_bin: Binary classification of the tourist traffic on the
  trail - high or low

- pseudotui: Continuous variable representing relative number of
  tourists per trail.

- value: Value 1, to be used for rasterization purposes

## Source

<https://register.geonorge.no/det-offentlige-kartgrunnlaget/n50-kartdata/ea192681-d039-42ec-b1bc-f3ce04c189ac>

## Examples

``` r
(f <- system.file("vector/hardanger_trails.gpkg", package = "oneimpact"))
#> [1] "/home/runner/work/_temp/Library/oneimpact/vector/hardanger_trails.gpkg"
v <- sf::st_read(f)
#> Reading layer `hardanger_trails' from data source 
#>   `/home/runner/work/_temp/Library/oneimpact/vector/hardanger_trails.gpkg' 
#>   using driver `GPKG'
#> Simple feature collection with 1444 features and 5 fields
#> Geometry type: MULTILINESTRING
#> Dimension:     XY
#> Bounding box:  xmin: 5168.794 ymin: 6613269 xmax: 194195.6 ymax: 6760287
#> Projected CRS: ETRS89 / UTM zone 33N
plot(v[1])

```
