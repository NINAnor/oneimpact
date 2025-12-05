# Create ZOI layer names as strings for data annotation

This function uses a vector of variable names and a vector of ZOI radii
to create a `data.frame` with all combinations of variable names and zoi
radii, to ease listing or accessing them from GRASS GIS or for setting
up a model formula.

## Usage

``` r
add_zoi_layers(layers, zoi_radius, name = NULL, pattern = NULL)
```

## Examples

``` r
# attaching zoi_radius to the end
name <- c("houses", "roads", "railways")
radii <- seq(100, 500, by = 100)
add_zoi_layers(name, radii)
#>    zoi_radius       layer
#> 1         100   houses100
#> 2         100    roads100
#> 3         100 railways100
#> 4         200   houses200
#> 5         200    roads200
#> 6         200 railways200
#> 7         300   houses300
#> 8         300    roads300
#> 9         300 railways300
#> 10        400   houses400
#> 11        400    roads400
#> 12        400 railways400
#> 13        500   houses500
#> 14        500    roads500
#> 15        500 railways500

# replacing a pattern by zoi_radius
themes <- c("houses_XXX", "private_cabins_XXX", "roads_XXX", "powerlines_XXX", "railways_XXX")
maps <- c("houses_XXX@my_mapset", "private_cabins_XXX@my_mapset", "roads_summer_XXX@my_mapset",
         "powerlines_XXX@my_mapset", "railway_XXX@my_mapset")

add_zoi_layers(layers = maps, zoi_radius = c(100, 500, 1000, 5000, 10000),
               name = themes, pattern = "XXX")
#>                    name zoi_radius                          layer
#> 1            houses_100        100           houses_100@my_mapset
#> 2    private_cabins_100        100   private_cabins_100@my_mapset
#> 3             roads_100        100     roads_summer_100@my_mapset
#> 4        powerlines_100        100       powerlines_100@my_mapset
#> 5          railways_100        100          railway_100@my_mapset
#> 6            houses_500        500           houses_500@my_mapset
#> 7    private_cabins_500        500   private_cabins_500@my_mapset
#> 8             roads_500        500     roads_summer_500@my_mapset
#> 9        powerlines_500        500       powerlines_500@my_mapset
#> 10         railways_500        500          railway_500@my_mapset
#> 11          houses_1000       1000          houses_1000@my_mapset
#> 12  private_cabins_1000       1000  private_cabins_1000@my_mapset
#> 13           roads_1000       1000    roads_summer_1000@my_mapset
#> 14      powerlines_1000       1000      powerlines_1000@my_mapset
#> 15        railways_1000       1000         railway_1000@my_mapset
#> 16          houses_5000       5000          houses_5000@my_mapset
#> 17  private_cabins_5000       5000  private_cabins_5000@my_mapset
#> 18           roads_5000       5000    roads_summer_5000@my_mapset
#> 19      powerlines_5000       5000      powerlines_5000@my_mapset
#> 20        railways_5000       5000         railway_5000@my_mapset
#> 21         houses_10000      10000         houses_10000@my_mapset
#> 22 private_cabins_10000      10000 private_cabins_10000@my_mapset
#> 23          roads_10000      10000   roads_summer_10000@my_mapset
#> 24     powerlines_10000      10000     powerlines_10000@my_mapset
#> 25       railways_10000      10000        railway_10000@my_mapset
```
