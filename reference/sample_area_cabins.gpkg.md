# Cabins vector data for the sample area

Dataset containing the location of tourist private cabins in Southern
Norway. It corresponds to some specific building types (object_type =
"Bygning", byggtyp_nbr = c("161", "162", "163")) form the public N50
dataset. The map was clipped for the sample area presented in the
`oneimpact` package.

## Format

A geopackage file. Projected CRS: ETRS89 / UTM zone 33N. The vector file
presents the following columns:

- cat: Line number, corresponding to the original dataset

- buildtype: Type of building (code) in the original dataset

- city: Code of the municipality where the cabin is located

- value: Value 1, to be used for rasterization purposes

## Source

<https://register.geonorge.no/det-offentlige-kartgrunnlaget/n50-kartdata/ea192681-d039-42ec-b1bc-f3ce04c189ac>

## See also

Maps for the sample area:  
Limits of sample area:
[sample_area.gpkg](https://ninanor.github.io/oneimpact/reference/sample_area.gpkg.md)  
Cabins:
[sample_area_cabins.tif](https://ninanor.github.io/oneimpact/reference/sample_area_cabins.tif.md),
[sample_area_cabins_count.tif](https://ninanor.github.io/oneimpact/reference/sample_area_cabins_count.tif.md)  
Roads:
[sample_area_roads.gpkg](https://ninanor.github.io/oneimpact/reference/sample_area_roads.gpkg.md),
[sample_area_roads.tif](https://ninanor.github.io/oneimpact/reference/sample_area_roads.tif.md)

## Examples

``` r
(f <- system.file("vector/sample_area_cabins.gpkg", package = "oneimpact"))
#> [1] "/home/runner/work/_temp/Library/oneimpact/vector/sample_area_cabins.gpkg"
sf::st_read(f)
#> Reading layer `sample_area_cabins' from data source 
#>   `/home/runner/work/_temp/Library/oneimpact/vector/sample_area_cabins.gpkg' 
#>   using driver `GPKG'
#> Simple feature collection with 6875 features and 4 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 146900.1 ymin: 6622822 xmax: 194694.6 ymax: 6658891
#> Projected CRS: ETRS89 / UTM zone 33N
# or
v <- terra::vect(f)
plot(v)

```
