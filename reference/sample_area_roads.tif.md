# Road raster data for the sample area

Raster data indicating pixels with public roads in Southern Norway.
Rasterized from the vector data from the Norwegian road dataset Elveg
1.0 with 100 m resolution and clipped for the study area presented in
the `oneimpact` package.

## Format

A Geotiff file. Projected CRS: ETRS89 / UTM zone 33N.

- 1: Presence of roads

- NA: No presence of roads

## Source

<https://kartkatalog.geonorge.no/metadata/elveg/ed1e6798-b3cf-48be-aee1-c0d3531da01a>

## See also

Maps for the sample area:  
Limits of sample area:
[sample_area.gpkg](https://ninanor.github.io/oneimpact/reference/sample_area.gpkg.md)  
Cabins:
[sample_area_cabins.gpkg](https://ninanor.github.io/oneimpact/reference/sample_area_cabins.gpkg.md),
[sample_area_cabins.tif](https://ninanor.github.io/oneimpact/reference/sample_area_cabins.tif.md),
[sample_area_cabins_count.tif](https://ninanor.github.io/oneimpact/reference/sample_area_cabins_count.tif.md)  
Roads:
[sample_area_roads.gpkg](https://ninanor.github.io/oneimpact/reference/sample_area_roads.gpkg.md)

## Examples

``` r
(f <- system.file("raster/sample_area_roads.tif", package = "oneimpact"))
#> [1] "/home/runner/work/_temp/Library/oneimpact/raster/sample_area_roads.tif"
r <- terra::rast(f)
plot(r)

```
