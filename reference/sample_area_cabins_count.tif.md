# Cabin count raster data for the sample area

Raster data indicating the number of tourist private cabins per pixel in
Norway. Cabins corresponds to some specific building types (object_type
= "Bygning", byggtyp_nbr = c("161", "162", "163")) form the public N50
dataset. The original data consisted of point vector data and were
rasterized with 100m resolution by counting the number of cabins in each
pixel. The raster was clipped for the study area presented in the
`oneimpact` package.

## Format

A Geotiff file. Projected CRS: ETRS89 / UTM zone 33N.

## Source

<https://register.geonorge.no/det-offentlige-kartgrunnlaget/n50-kartdata/ea192681-d039-42ec-b1bc-f3ce04c189ac>

## See also

Maps for the sample area:  
Limits of sample area:
[sample_area.gpkg](https://ninanor.github.io/oneimpact/reference/sample_area.gpkg.md)  
Cabins:
[sample_area_cabins.gpkg](https://ninanor.github.io/oneimpact/reference/sample_area_cabins.gpkg.md),
[sample_area_cabins.tif](https://ninanor.github.io/oneimpact/reference/sample_area_cabins.tif.md)  
Roads:
[sample_area_roads.gpkg](https://ninanor.github.io/oneimpact/reference/sample_area_roads.gpkg.md),
[sample_area_roads.tif](https://ninanor.github.io/oneimpact/reference/sample_area_roads.tif.md)

## Examples

``` r
(f <- system.file("raster/sample_area_cabins_count.tif", package = "oneimpact"))
#> [1] "/home/runner/work/_temp/Library/oneimpact/raster/sample_area_cabins_count.tif"
r <- terra::rast(f)
plot(r)

```
