# libraries
library(rgrass)
library(terra)

# Load raster data
f <- system.file("raster/sample_area_cabins.tif", package = "oneimpact")
cabins <- terra::rast(f)

# connect to grass gis and create grass location
# For linux or within OSGeo4W shell
grassdir <- system("grass --config path", intern = TRUE)
# grassdir <- system("grass78 --config path", intern = TRUE) # for GRASS 7.8
# If you used the standalone installer in Windows
# grassdir <- "C:\\Programs\\GRASS GIS 7.8" # Correct if the path GRASS version or path is different

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
rgrass::write_RAST(cabins, "cabins", flags = "o")

# binarize the input map

# map with only 1
cabins_bin1_name <- grass_binarize("cabins", output = "cabins_bin1",
                                   breaks = 1, overwrite = T)
# map with 0, 1
cabins_bin2_name <- grass_binarize("cabins", output = "cabins_bin2",
                                   breaks = 1, null = 0, overwrite = T)

# visualize
cabins_bin1_2 <- rgrass::read_RAST(c(cabins_bin1_name, cabins_bin2_name),
                                   return_format = "terra", NODATA = 255)
plot(cabins_bin1_2, main = c("Binarized map keeping null", "Binarized map setting null to 0"))

#-------
# binarize the map with multiple break values

# first create a continuous map
cont_map_name <- calc_zoi_nearest("cabins_bin1", radius = 1000,
                                  type = "exp_decay",
                                  where = "GRASS", overwrite = TRUE)
# binarize
cabins_bin2vals_name <- grass_binarize(cont_map_name, output = "cabins_zoi1000_bin",
                                       breaks = c(0.3, 0.8), overwrite = T)
# visualize
cabins_bin2vals <- rgrass::read_RAST(c(cont_map_name, cabins_bin2vals_name),
                                     return_format = "terra", NODATA = 255)
plot(cabins_bin2vals,
     main = c("Original map",
              "Binarized map, break = 0.3",
              "Binarized map, break = 0.8"))
