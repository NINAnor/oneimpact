# Running calc_zoi_nearest through GRASS GIS
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

# define map name within GRASS GIS
cabins_g <- "cabins_example"
# add file to GRASS GIS mapset
rgrass::write_RAST(cabins, cabins_g, flags = c("overwrite", "o"))

# check
terra::plot(cabins, col = "black",
            main = "Map of tourist cabins")

#---
# define region in GRASS GIS
rgrass::execGRASS("g.region", raster = cabins_g,
                  flags = "p")

# Input map name within GRASS GIS
cabins_g

# Exponential decay ZoI=1000m
expdecay_name <- calc_zoi_nearest(cabins_g, type = "exp_decay",
                                  radius = 1000,
                                  where = "GRASS",
                                  g_verbose = FALSE, g_overwrite = TRUE)

# Bartlett decay ZoI=1000m
bartlett_name <- calc_zoi_nearest(cabins_g, type = "bartlett",
                                  radius = 1000,
                                  where = "GRASS", g_verbose = FALSE, g_overwrite = TRUE)

# Threshold influence ZoI = 1000m
threshold_name <- calc_zoi_nearest(cabins_g, type = "threshold",
                                   radius = 1000,
                                   where = "GRASS", g_verbose = FALSE, g_overwrite = TRUE)

# Gaussian influence ZoI = 1000m
gaussian_name <- calc_zoi_nearest(cabins_g, type = "Gauss",
                                  radius = 1000,
                                  where = "GRASS", g_verbose = FALSE, g_overwrite = TRUE)

# Log-distance
log_name <- calc_zoi_nearest(cabins_g, type = "log", log_base = 10,
                             where = "GRASS",
                             g_verbose = FALSE, g_overwrite = TRUE)

# Euclidean
euclidean_name <- calc_zoi_nearest(cabins_g, type = "euclidean",
                                   where = "GRASS",
                                   g_verbose = FALSE, g_overwrite = TRUE)

(all_names <- c(euclidean_name, log_name, expdecay_name,
                bartlett_name, threshold_name, gaussian_name))

# visualize
cabins_zoi_nearest <- rgrass::read_RAST(all_names, return_format = "terra")

title_plot <- c("Euclidean distance", "Log distance (base 10)",
                "Exponential ZoI 1000m", "Bartlett ZoI 1000m",
                "Threshold ZoI 1000m", "Gaussian ZoI 1000m")
terra::plot(cabins_zoi_nearest, main = title_plot)

# remove rasters created
# to_remove_rast <- c(all_names)
# rgrass::execGRASS("g.remove", type = "rast", name = to_remove_rast, flags = "f")
