# Running calc_influence_nearest through GRASS GIS
library(rgrass7)
library(raster)
library(terra)
library(sp)
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
terra::plot(cabins, col = "black")

#---
# define region in GRASS GIS
rgrass7::execGRASS("g.region", raster = cabins_g,
                   flags = "p")

# Input map name within GRASS GIS
cabins_g

# Euclidean
euclidean_name <- calc_influence_nearest(cabins_g, where = "GRASS",
                                         quiet = T, overwrite = T)
# Log
log_name <- calc_influence_nearest(cabins_g, type = "log", log_base = 10,
                                   where = "GRASS", quiet = T, overwrite = T)
# Exponential decay ZoI=1000m
expdecay_name <- calc_influence_nearest(cabins_g, type = "exp_decay", zoi = 1000,
                                        where = "GRASS", quiet = T, overwrite = T)
# Bartlett decay ZoI=1000m
bartlett_name <- calc_influence_nearest(cabins_g, type = "bartlett", zoi = 1000,
                                        where = "GRASS", quiet = T, overwrite = T)
# Threshold influence ZoI = 1000m
threshold_name <- calc_influence_nearest(cabins_g, type = "threshold", zoi = 1000,
                                         where = "GRASS", quiet = T, overwrite = T)
# Gaussian influence ZoI = 1000m
gaussian_name <- calc_influence_nearest(cabins_g, type = "Gauss", zoi = 1000,
                                         where = "GRASS", quiet = T, overwrite = T)


(all_names <- c(euclidean_name, log_name, expdecay_name, bartlett_name, threshold_name, gaussian_name))

# visualize
cabins_influence_nearest <- readRAST(all_names) %>%
  raster::stack() %>%
  terra::rast()

title_plot <- c("Euclidean distance", "Log distance (base 10)",
                "Exponential decay 1000m", "Bartlett decay 1000m",
                "Threshold influence 1000m", "Gaussian influence 1000m")
terra::plot(cabins_influence_nearest, main = title_plot)

# remove rasters created
to_remove_vect <- c(test_region_name, cabins_vect_name)
to_remove_rast <- c(all_names)
# rgrass7::execGRASS("g.remove", type = "vect", name = to_remove_vect, flags = "f")
# rgrass7::execGRASS("g.remove", type = "rast", name = to_remove_rast, flags = "f")
