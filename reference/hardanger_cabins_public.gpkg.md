# Public cabins vector data for the Hardangervidda wild reindeer area in Norway

Dataset containing the location of large, public DNT cabins and mountain
hotels in the surroundings of the Hardangervidda wild reindeer area.
Retrieved from the public N50 dataset.

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
v <- terra::vect(f)
plot(v)

```
