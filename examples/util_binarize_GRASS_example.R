# libraries
library(rgrass7)
library(raster)
library(terra)
library(dplyr)

# connect to grass gis 7.8
grassdir <- system("grass78 --config path", intern = T)
gisDB <- "/data/grass"
loc <- "ETRS_33N/"
ms <- "u_bb_cuminf"
rgrass7::initGRASS(gisBase = grassdir,
                   home = tempdir(),
                   override = T,
                   gisDbase = gisDB,
                   location = loc,
                   mapset = ms)

# input map (not binarized)
rgrass7::use_sp()
name_var <- "private_cabins_sub"

# binarize the input map

# map with only 1
cabins_bin1_name <- util_binarize_GRASS(name_var, output = "cabins_bin1",
                                        breaks = 1, overwrite = T)
# map with 0, 1
cabins_bin2_name <- util_binarize_GRASS(name_var, output = "cabins_bin2",
                                        breaks = 1, null = 0, overwrite = T)

# visualize
cabins_bin1_2 <- rgrass7::readRAST(c(cabins_bin2_name, cabins_bin1_name)) %>%
  raster::stack() %>%
  terra::rast()
plot(cabins_bin1_2, main = c("Binarized map setting null to 0", "Binarized map keeping null"))

# binarize the map with multiple break values

# first create a continuous map
cont_map_name <- calc_influence_nearest(name_var, zoi = 1000, transform = "exp_decay",
                                        where = "GRASS", overwrite = T)
# binarize
cabins_bin2vals_name <- util_binarize_GRASS(cont_map_name, output = "cabins_bin",
                                            breaks = c(0.3, 0.5), overwrite = T)
# visualize
cabins_bin2vals <- rgrass7::readRAST(c(cont_map_name, cabins_bin2vals_name)) %>%
  raster::stack() %>%
  terra::rast()
plot(cabins_bin2vals,
     main = c("Original map",
              "Binarized map, break = 0.3",
              "Binarized map, break = 0.5"))
