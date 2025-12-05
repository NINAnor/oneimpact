# Rasterizes a vector counting the number of features in each pixel

This function rasterizes a vector file in GRASS GIS, counting the number
of vector features within each pixel. The function uses the resolution
and extent already set as the GRASS GIS mapset's computational region.
If `input_as_region` is set to `TRUE`, the extent of the vector map `x`
is used reset as the region. The resolution, though, continues the same
set up earlier through `g.region`.

## Usage

``` r
grass_v2rast_count(
  x,
  output = paste0(x, "_count"),
  column = NULL,
  input_as_region = FALSE,
  align = NULL,
  remove_intermediate = TRUE,
  verbose = FALSE,
  overwrite = FALSE,
  ...
)
```

## Arguments

- x:

  Input vector map.

- output:

  Output map name.

- column:

  `[chracter(1)=NULL`  
  Default is `NULL`. If not `NULL`, the name of a column in the input
  vector `x` that corresponds to the column to be summed to count the
  number of features in each pixel of the output raster map. If `NULL`,
  this column is created in a temporary vector, with all values equal 1.

- input_as_region:

  `[logical(1)=FALSE]`  
  Default is FALSE. Whether the GRASS GIS computational region should be
  set within the function (to the extent of `x`) or not. If `FALSE`, the
  current computational region is used.

- align:

  `[character(1)=NULL]`  
  Name of a raster map with which to align the computational region to
  produce the output map.

- verbose:

  `[logical(1)=FALSE]`  
  Should messages of the computation steps be printed in the prompt
  along the computation?

- overwrite:

  `[logical(1)]`  
  Whether the output maps should be overwriten (flag
  `overwrite = TRUE`).

## Value

A raster map with the count of features within each pixel. The map is
written within the GRASS GIS mapset. In R, the output is a string with
the name of this map.

## Examples

``` r
# not run
if (FALSE) { # \dontrun{

  # libraries
  library(rgrass)
  library(terra)

  # Load vector data
  f <- system.file("vector/sample_area_cabins.gpkg", package = "oneimpact")
  cabins_vect <- terra::vect(f)
  # Load referemce raster data - for creation of GRASS project
  f <- system.file("raster/sample_area_cabins.tif", package = "oneimpact")
  cabins <- terra::rast(f)

  # connect to grass gis and create grass location
  # For linux or within OSGeo4W shell
  grassdir <- system("grass --config path", intern = TRUE)
  # grassdir <- system("grass78 --config path", intern = TRUE) # for GRASS 7.8
  # If you used the standalone installer in Windows
  # grassdir <- "C:\Programs\GRASS GIS 7.8" # Correct if the path GRASS version or path is different

  gisDB <- "." # create location and mapset in the working directory
  loc <- "ETRS_33N/" # name of the location
  ms <- "PERMANENT" # name of the mapset
  rgrass::initGRASS(gisBase = grassdir,
                    SG = cabins, # use map to define location projection
                    home = tempdir(),
                    override = TRUE,
                    gisDbase = gisDB,
                    location = loc,
                    mapset = ms)

  # add map to GRASS
  rgrass::write_VECT(cabins_vect, vname = "cabins_vect", flags = "o")

  # set region
  rgrass::execGRASS("g.region", vector = cabins_vect, res = "100", flags = c("p"))

  # rasterize with count, creating a new temp_vector
  cabins_count_name <- grass_v2rast_count("cabins_vect", output = "cabins_count",
                                          verbose = TRUE, overwrite = TRUE)

  # rasterize with count, without creating a temporary vector
  cabins_count_name <- grass_v2rast_count(cabins_vect_name, output = "cabins_count",
                                          column = "value",
                                          verbose = TRUE, overwrite = TRUE)

  # visualize
  rgrass::read_RAST(cabins_count_name, return_format = "terra")
  plot(main = "Number of private cabins")

  # remove rasters created
  # to_remove_vect <- c(test_region_name, cabins_vect_name)
  # to_remove_rast <- c(cabins_count_name)
  # rgrass::execGRASS("g.remove", type = "vect", name = to_remove_vect, flags = "f")
  # rgrass::execGRASS("g.remove", type = "rast", name = to_remove_rast, flags = "f")

} # }
```
