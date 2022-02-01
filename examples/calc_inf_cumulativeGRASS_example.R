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

# define region
##########################
# change it here when these data are within the R package

# I defined a region_test_influence through the GRASS GUI and using v.in.region output=region_test_influence
rgrass7::execGRASS("g.region", parameters = list(vector = "region_test_influence", align = "private_cabins_rast"),
                   flags = "p")

# show input map and region here
rgrass7::use_sp()
cabins <- rgrass7::readRAST(c("private_cabins_rast_sub")) %>%
  raster::raster() %>%
  terra::rast()

terra::plot(cabins, col = "black")

rgrass7::use_sf()
study_area <- rgrass7::readVECT(c("region_test_influence"))


# Input map name within GRASS GIS
name_var <- "private_cabins_sub"

# exp_decay R
expR <- calc_influence_cumulative(x = cabins, zoi = 1000, type = "exp_decay")
plot(expR)

expR2 <- calc_influence_cumulative(x = cabins, zoi = c(1000, 3000), type = "exp_decay", normalize = F)
plot(expR2)

# normalization
# no: 0 to 100 cabins
# in the end: sum = 1
# before the end: where we have max number of cabins = 1

expG_name <- calc_influence_cumulative_GRASS(x = name_var, zoi = 1000, type = "exp_decay",
                                             overwrite = T, quiet = F)

expG_name2 <- calc_influence_cumulative_GRASS(x = name_var, zoi = c(1000, 2000), type = "exp_decay",
                                             overwrite = T, quiet = F)
# visualize
expG <- rgrass7::readRAST(expG_name2) %>%
  raster::raster() %>%
  terra::rast()

plot(expG)
plot(c(expR, expG))

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
