# not run
\dontrun{

# libraries
library(rgrass7)
library(terra)

# Load vector data
f <- system.file("vector/sample_area_cabins.gpkg", package = "oneimpact")
cabins_vect <- terra::vect(f)
# Load referemce raster data - for creation of GRASS project
f <- system.file("raster/sample_area_cabins.tif", package = "oneimpact")
cabins <- terra::rast(f)

# connect to grass gis 7.8 and create grass location
# For linux or within OSGeo4W shell
grassdir <- system("grass78 --config path", intern = TRUE)
# If you used the standalone installer in Windows
# grassdir <- "C:\\Programs\\GRASS GIS 7.8" # Correct if the path is different

gisDB <- "." # create location and mapset in the working directory
loc <- "ETRS_33N/" # name of the location
ms <- "PERMANENT" # name of the mapset
rgrass7::initGRASS(gisBase = grassdir,
                   SG = cabins, # use map to define location projection
                   home = tempdir(),
                   override = TRUE,
                   gisDbase = gisDB,
                   location = loc,
                   mapset = ms)

# add map to GRASS
rgrass7::write_VECT(cabins_vect, vname = "cabins_vect", flags = "o")

# set region
rgrass7::execGRASS("g.region", vector = cabins_vect, res = "100", flags = c("p"))

# rasterize with count, creating a new temp_vector
cabins_count_name <- grass_v2rast_count("cabins_vect", output = "cabins_count",
                                        verbose = TRUE, overwrite = TRUE)

# rasterize with count, without creating a temporary vector
cabins_count_name <- grass_v2rast_count(cabins_vect_name, output = "cabins_count",
                                        column = "value",
                                        verbose = TRUE, overwrite = TRUE)

# visualize
rgrass7::read_RAST(cabins_count_name, return_format = "terra")
plot(main = "Number of private cabins")

# remove rasters created
# to_remove_vect <- c(test_region_name, cabins_vect_name)
# to_remove_rast <- c(cabins_count_name)
# rgrass7::execGRASS("g.remove", type = "vect", name = to_remove_vect, flags = "f")
# rgrass7::execGRASS("g.remove", type = "rast", name = to_remove_rast, flags = "f")

}
