# Sample area: a polygon vector data

Dataset containing the limits of an arbitrary study area in Southern
Norway, used for illustrative purposes.

## Format

A geopackage file. Projected CRS: ETRS89 / UTM zone 33N.

## See also

Maps for the sample area:  
Cabins:
[sample_area_cabins.gpkg](https://ninanor.github.io/oneimpact/reference/sample_area_cabins.gpkg.md),
[sample_area_cabins.tif](https://ninanor.github.io/oneimpact/reference/sample_area_cabins.tif.md),
[sample_area_cabins_count.tif](https://ninanor.github.io/oneimpact/reference/sample_area_cabins_count.tif.md)  
Roads:
[sample_area_roads.gpkg](https://ninanor.github.io/oneimpact/reference/sample_area_roads.gpkg.md),
[sample_area_roads.tif](https://ninanor.github.io/oneimpact/reference/sample_area_roads.tif.md)

## Examples

``` r
(f <- system.file("vector/sample_area.gpkg", package = "oneimpact"))
#> [1] "/home/runner/work/_temp/Library/oneimpact/vector/sample_area.gpkg"
sf::st_read(f)
#> Reading layer `study_area' from data source 
#>   `/home/runner/work/_temp/Library/oneimpact/vector/sample_area.gpkg' 
#>   using driver `GPKG'
#> Simple feature collection with 1 feature and 1 field
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 146900 ymin: 6622800 xmax: 194700 ymax: 6658900
#> Projected CRS: ETRS89 / UTM zone 33N
# or
v <- terra::vect(f)
plot(v)

```
