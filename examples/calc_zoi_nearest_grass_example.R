# Running calc_zoi_nearest through GRASS GIS
library(rgrass7)
library(terra)
library(sp)
library(dplyr)

# Load raster data
f <- system.file("raster/cabins.tif", package = "oneimpact")
cabins <- terra::rast(f)

# connect to grass gis 7.8 and create grass location
grassdir <- system("grass78 --config path", intern = T)
gisDB <- "." # create location and mapset in the working directory
loc <- "ETRS_33N/" # name of the location
ms <- "PERMANENT" # name of the mapset
rgrass7::initGRASS(gisBase = grassdir,
                   SG = cabins, # use map to define location projection
                   home = tempdir(),
                   override = T,
                   gisDbase = gisDB,
                   location = loc,
                   mapset = ms)


# define map name within GRASS GIS
cabins_g <- "cabins_example"
# add file to GRASS GIS mapset
rgrass7::write_RAST(cabins, cabins_g, flags = "overwrite")

# check
terra::plot(cabins, col = "black",
            main = "Map of tourist cabins")

#---
# define region in GRASS GIS
rgrass7::execGRASS("g.region", raster = cabins_g,
                   flags = "p")

# Input map name within GRASS GIS
cabins_g

# Euclidean
euclidean_name <- calc_zoi_nearest(cabins_g, type = "euclidean",
                                   where = "GRASS",
                                   quiet = T, overwrite = T)
# Log
log_name <- calc_zoi_nearest(cabins_g, type = "log", log_base = 10,
                             where = "GRASS", quiet = T, overwrite = T)
# Exponential decay ZoI=1000m
expdecay_name <- calc_zoi_nearest(cabins_g, type = "exp_decay",
                                  radius = 1000,
                                  where = "GRASS", quiet = T, overwrite = T)
# Bartlett decay ZoI=1000m
bartlett_name <- calc_zoi_nearest(cabins_g, type = "bartlett",
                                  radius = 1000,
                                  where = "GRASS", quiet = T, overwrite = T)

# Threshold influence ZoI = 1000m
threshold_name <- calc_zoi_nearest(cabins_g, type = "threshold",
                                   radius = 1000,
                                   where = "GRASS", quiet = T, overwrite = T)

# Gaussian influence ZoI = 1000m
gaussian_name <- calc_zoi_nearest(cabins_g, type = "Gauss",
                                  radius = 1000,
                                  where = "GRASS", quiet = T, overwrite = T)


(all_names <- c(euclidean_name, log_name, expdecay_name,
                bartlett_name, threshold_name, gaussian_name))

# visualize
cabins_zoi_nearest <- rgrass7::read_RAST(all_names, return_format = "terra")

title_plot <- c("Euclidean distance", "Log distance (base 10)",
                "Exponential ZoI 1000m", "Bartlett ZoI 1000m",
                "Threshold ZoI 1000m", "Gaussian ZoI 1000m")
terra::plot(cabins_zoi_nearest, main = title_plot)

# remove rasters created
to_remove_vect <- c(test_region_name, cabins_vect_name)
to_remove_rast <- c(all_names)
# rgrass7::execGRASS("g.remove", type = "vect", name = to_remove_vect, flags = "f")
# rgrass7::execGRASS("g.remove", type = "rast", name = to_remove_rast, flags = "f")
