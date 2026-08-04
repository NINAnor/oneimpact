# Public roads vector data for the reindeer area

Dataset containing the location of large, public roads in Southern
Norway, within the reindeer management area of Setesdal Austhei.
Retrieved from the Norwegian road dataset Elveg 1.0.

## Format

A geopackage file. Projected CRS: ETRS89 / UTM zone 33N. The vector file
presents the following columns:

- id: Line number, corresponding to the original dataset

- name: Local name of the road

- publ_priv: Whether the road is public or private

- traffic_bin: Binary classification of the traffic on the road - high
  or low

- name_area: Name of the reindeer management area where the road is
  located

- traffic_bin: Value 1, to be used for rasterization purposes

## Source

<https://kartkatalog.geonorge.no/metadata/elveg-20/77944f7e-3d75-4f6d-ae04-c528cc72e8f6>

## Examples

``` r
(f <- system.file("vector/reindeer_roads_public.gpkg", package = "oneimpact"))
#> [1] "/home/runner/work/_temp/Library/oneimpact/vector/reindeer_roads_public.gpkg"
v <- sf::st_read(f)
#> Reading layer `reindeer_roads_public' from data source 
#>   `/home/runner/work/_temp/Library/oneimpact/vector/reindeer_roads_public.gpkg' 
#>   using driver `GPKG'
#> Simple feature collection with 2862 features and 6 fields
#> Geometry type: MULTILINESTRING
#> Dimension:     XY
#> Bounding box:  xmin: 63605.12 ymin: 6513202 xmax: 130683.5 ymax: 6645648
#> Projected CRS: ETRS89 / UTM zone 33N
plot(v[1])

```
