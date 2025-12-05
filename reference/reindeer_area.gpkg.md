# Reindeer area: a polygon vector data for the Setesdal Austhei reindeer herding area

Dataset containing the limits of the reindeer management area of
Setesdal Austhei in Southern Norway, as defined in Panzacchi et al.
(2015). Note this area is slightly larger than the boundaries used for
reindeer management in Norway.

## Format

A geopackage file. Projected CRS: ETRS89 / UTM zone 33N.

## Source

Panzacchi, M., Van Moorter, B., Strand, O., Loe, L. E., & Reimers, E.
(2015). Searching for the fundamental niche using individual-based
habitat selection modelling across populations. Ecography, 38(7),
659–669. <https://doi.org/10.1111/ecog.01075>

## Examples

``` r
(f <- system.file("vector/reindeer_area.gpkg", package = "oneimpact"))
#> [1] "/home/runner/work/_temp/Library/oneimpact/vector/reindeer_area.gpkg"
sf::st_read(f)
#> Reading layer `reindeer_area' from data source 
#>   `/home/runner/work/_temp/Library/oneimpact/vector/reindeer_area.gpkg' 
#>   using driver `GPKG'
#> Simple feature collection with 1 feature and 3 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 63586.78 ymin: 6513160 xmax: 130692.6 ymax: 6645648
#> Projected CRS: ETRS89 / UTM zone 33N
# or
v <- terra::vect(f)
plot(v)

```
