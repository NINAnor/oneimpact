# Polygon delimiting the study area around the Hardangervidda wild reindeer area

Dataset containing the polygon of the Hardangervidda area. It is
different from the official wild reindeer area and is delimited by the
main roads and barriers in the surroundings of the Hardangervidda
plateau in Southern Norway,

## Format

A geopackage file. Projected CRS: ETRS89 / UTM zone 33N. The vector file
presents the following columns:

- reindeer_areas_id: ID of the area

- name_area: Name of the wild reindeer area

## Source

Niebuhr, B. B., Van Moorter, B., Stien, A., Tveraa, T., Strand, O.,
Langeland, K., Sandström, P., Alam, M., Skarin, A., & Panzacchi, M.
(2023). Estimating the cumulative impact and zone of influence of
anthropogenic features on biodiversity. Methods in Ecology and
Evolution. https://doi.org/10.1111/2041-210X.14133

## Examples

``` r
(f <- system.file("vector/hardanger_polygon.gpkg", package = "oneimpact"))
#> [1] "/home/runner/work/_temp/Library/oneimpact/vector/hardanger_polygon.gpkg"
v <- terra::vect(f)
plot(v)

```
