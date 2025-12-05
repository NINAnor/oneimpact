# Cabins vector data for the reindeer area

Dataset containing the location of tourist private cabins is Southern
Norway, within the reindeer management area of Setesdal Austhei. It
corresponds to some specific building types (object_type = "Bygning",
byggtyp_nbr = c("161", "162", "163")) from the public N50 dataset.

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
(f <- system.file("vector/reindeer_cabins.gpkg", package = "oneimpact"))
#> [1] "/home/runner/work/_temp/Library/oneimpact/vector/reindeer_cabins.gpkg"
sf::st_read(f)
#> Reading layer `reindeer_cabins' from data source 
#>   `/home/runner/work/_temp/Library/oneimpact/vector/reindeer_cabins.gpkg' 
#>   using driver `GPKG'
#> Simple feature collection with 5605 features and 4 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 63646.4 ymin: 6513651 xmax: 130352.7 ymax: 6644403
#> Projected CRS: ETRS89 / UTM zone 33N
# or
v <- terra::vect(f)
plot(v)

```
