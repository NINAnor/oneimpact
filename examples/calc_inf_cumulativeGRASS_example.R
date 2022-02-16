# Running calc_influence_cumulative through GRASS GIS
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

# Load raster data
f <- system.file("raster/cabins.tif", package = "oneimpact")
cabins_sp <- raster::raster(f) %>%
  as("SpatialPixelsDataFrame")
# define map name within GRASS GIS
cabins_g <- "private_cabins_sub"
# add file to GRASS GIS mapset
rgrass7::use_sp()
rgrass7::writeRAST(cabins_sp, cabins_g, overwrite = TRUE)

# check
cabins <- cabins_sp %>%
  raster::raster() %>%
  terra::rast()
terra::plot(cabins, col = "black",
            main = "Map of cabins")

#---
# define region in GRASS GIS
rgrass7::execGRASS("g.region", raster = cabins_g,
                   flags = "p")

# Input map name within GRASS GIS - binary map
cabins_bin_g <- util_binarize_GRASS(cabins_g, output = "private_cabins_bin",
                                    null = 0, overwrite = TRUE)

# check input
cabins_bin <- rgrass7::readRAST(cabins_bin_g) %>%
  raster::raster() %>%
  terra::rast()
plot(cabins_bin, col = c("lightyellow", "black"),
     main = "Binarized map of cabins")

# Exponential decay
exp_name <- calc_influence_cumulative(x = cabins_bin_g, zoi = 1000,
                                      zoi_decay_threshold = 0.01, type = "exp_decay",
                                      where = "GRASS", overwrite = T, quiet = F)
# Bartlett decay
barlett_name <- calc_influence_cumulative(x = cabins_bin_g, zoi = 1000, type = "bartlett",
                                      where = "GRASS", overwrite = T, quiet = F)
# Gaussian decay
gauss_name <- calc_influence_cumulative(x = cabins_bin_g, zoi = 1000,
                                        zoi_decay_threshold = 0.01, type = "Gauss",
                                        where = "GRASS", overwrite = T, quiet = F)
# Threshold decay (circle, step)
threshold_name <- calc_influence_cumulative(x = cabins_bin_g, zoi = 1000, type = "threshold",
                                            where = "GRASS", overwrite = T, quiet = F)

(all_names <- c(exp_name, barlett_name, gauss_name, threshold_name))

# visualize
cabins_influence_cumulative <- rgrass7::readRAST(all_names) %>%
  raster::stack() %>%
  terra::rast()

title_plot <- c("Exponential decay 1000m", "Bartlett decay 1000m",
                "Gaussian decay 1000m", "Threshold decay 1000m")
terra::plot(cabins_influence_cumulative, main = title_plot)

# remove rasters created
to_remove_vect <- c(test_region_name, cabins_vect_name)
to_remove_rast <- c(all_names)
# rgrass7::execGRASS("g.remove", type = "vect", name = to_remove_vect, flags = "f")
# rgrass7::execGRASS("g.remove", type = "rast", name = to_remove_rast, flags = "f")
