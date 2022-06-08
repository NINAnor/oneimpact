# Running calc_zoi_cumulative through GRASS GIS
library(rgrass7)
library(terra)
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
rgrass7::write_RAST(cabins, cabins_g, flags = c("o", "overwrite"))

# check
terra::plot(cabins, col = "black",
            main = "Map of tourist cabins")

#---
# define region in GRASS GIS
rgrass7::execGRASS("g.region", raster = cabins_g,
                   flags = "p")

# Input map name within GRASS GIS - binary map
cabins_bin_g <- util_binarize_grass(cabins_g, breaks = 1, output = "cabins_example_bin",
                                    null = 0, overwrite = TRUE)

# check input
cabins_bin <- rgrass7::read_RAST("cabins_example_bin", return_format = "terra", NODATA = 255)

plot(cabins_bin, col = c("lightyellow", "black"),
     main = "Binarized map of cabins")

# Exponential decay
exp_name <- calc_zoi_cumulative(x = cabins_bin_g,
                                zoi_radius = 1000, zoi_limit = 0.01,
                                type = "exp_decay",
                                where = "GRASS",
                                overwrite = T, quiet = F)
# Bartlett decay
barlett_name <- calc_zoi_cumulative(x = cabins_bin_g, zoi_radius = 1000,
                                    type = "bartlett",
                                    where = "GRASS",
                                    overwrite = T, quiet = F)
# Gaussian decay
gauss_name <- calc_zoi_cumulative(x = cabins_bin_g,
                                  zoi_radius = 1000, zoi_limit = 0.01,
                                  type = "Gauss",
                                  where = "GRASS",
                                  overwrite = T, quiet = F)

# Threshold decay (circle, step)
threshold_name <- calc_zoi_cumulative(x = cabins_bin_g, zoi_radius = 1000,
                                      type = "threshold",
                                      where = "GRASS",
                                      overwrite = T, quiet = F)

(all_names <- c(exp_name, barlett_name, gauss_name, threshold_name))

# visualize
cabins_zoi_cumulative <- rgrass7::read_RAST(all_names, return_format = "terra")

title_plot <- c("Exponential decay 1000m", "Bartlett decay 1000m",
                "Gaussian decay 1000m", "Threshold decay 1000m")
terra::plot(cabins_zoi_cumulative, main = title_plot)

#---
# calculate density vs cumulative ZoI
exp_name_d <- calc_zoi_cumulative(x = cabins_bin_g,
                                  zoi_radius = 1000, zoi_limit = 0.01,
                                  type = "exp_decay", output_type = "density",
                                  where = "GRASS",
                                  overwrite = T, quiet = F)

cabins_density <- rgrass7::read_RAST(exp_name_d, return_format = "terra")

terra::plot(c(cabins_zoi_cumulative[[1]], cabins_density),
            main = c("Cumulative ZoI", "Density"))

#---
# Using 'r.resamp.filter' algorithm

# rectangle
rectangle_resamp_filt <- calc_zoi_cumulative(x = cabins_bin_g,
                                             zoi_radius = 1000,
                                             type = "box",
                                             output_type = "density",
                                             where = "GRASS",
                                             module = "r.resamp.filter",
                                             overwrite = T, quiet = F)
rgrass7::read_RAST(rectangle_resamp_filt, return_format = "terra") %>%
  plot(main = "Rectangle ZoI 1000m")

# bartlett
bartlett_resamp_filt <- calc_zoi_cumulative(x = cabins_bin_g,
                                            zoi_radius = 1000,
                                            type = "bartlett",
                                            output_type = "cumulative_zoi",
                                            where = "GRASS",
                                            module = "r.resamp.filter",
                                            overwrite = T, quiet = F)
rgrass7::read_RAST(bartlett_resamp_filt, return_format = "terra") %>%
  plot(main = "Bartlett ZoI 1000m")

# Gaussian - to be implemented!
if(FALSE) {
  gauss_resamp_filt <- calc_zoi_cumulative(x = cabins_bin_g,
                                           zoi_radius = "1000,3000",
                                           type = "gauss,box",
                                           output_type = "cumulative_zoi",
                                           where = "GRASS",
                                           module = "r.resamp.filter",
                                           overwrite = T, quiet = F)
  rgrass7::read_RAST(bartlett_resamp_filt, return_format = "terra") %>%
    plot()
}

# remove rasters created
to_remove_rast <- unique(c(all_names, exp_name_d,
                           rectangle_resamp_filt, bartlett_resamp_filt))
# rgrass7::execGRASS("g.remove", type = "vect", name = to_remove_vect, flags = "f")
# rgrass7::execGRASS("g.remove", type = "rast", name = to_remove_rast, flags = "f")
