# Running calc_influence_nearest through GRASS GIS
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

# define region
##########################
# change it here when these data are within the R package

# I defined a region_test_influence through the GRASS GUI and using v.in.region output=region_test_influence
rgrass7::execGRASS("g.region", parameters = list(vector = "region_test_influence", align = "private_cabins_rast"),
                   flags = "p")

# show input map and region here
rgrass7::use_sp()
cabins <- readRAST(c("private_cabins_rast_sub")) %>% 
  raster::raster() %>% 
  terra::rast()

terra::plot(cabins, col = "black")

# Input map name within GRASS GIS
name_var <- "private_cabins_sub"
# Euclidean
euclidean_name <- calc_influence_nearest(name_var, where = "GRASS", 
                                         quiet = T, overwrite = T)
# Log
log_name <- calc_influence_nearest(name_var, transform = "log", log_base = 10,
                                   where = "GRASS", quiet = T, overwrite = T)
# Exponential decay ZoI=1000m
expdecay_name <- calc_influence_nearest(name_var, transform = "exp_decay", zoi = 1000,
                                        where = "GRASS", quiet = T, overwrite = T)
# Bartlett decay ZoI=1000m
bartlett_name <- calc_influence_nearest(name_var, transform = "bartlett", zoi = 1000,
                                        where = "GRASS", quiet = T, overwrite = T)

(all_names <- c(euclidean_name, log_name, expdecay_name, bartlett_name))

# visualize
cabins_influence_nearest <- readRAST(all_names) %>% 
  raster::stack() %>% 
  terra::rast()

title_plot <- c("Euclidean distance", "Log distance (base 10)", 
                "Exponential decay 1000m", "Bartlett decay 1000m")
terra::plot(cabins_influence_nearest, main = title_plot)
