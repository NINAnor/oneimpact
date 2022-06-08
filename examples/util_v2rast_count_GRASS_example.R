# libraries
library(rgrass7)
library(raster)
library(terra)
library(dplyr)
library(sf)

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

# input region
rgrass7::use_sf()
(s <- system.file("vector/study_area.gpkg", package = "oneimpact"))
test_region_name <- "region_test_oneimpact_pkg"
sf::st_read(s) %>%
  rgrass7::writeVECT(vname = test_region_name, v.in.ogr_flags = "overwrite")

# set region
rgrass7::execGRASS("g.region", vector = test_region_name, res = "100", flags = c("p"))

# input vector
(c <- system.file("vector/cabins_vect.gpkg", package = "oneimpact"))
cabins_vect_name <- "private_cabins_vect"
sf::st_read(c) %>%
  rgrass7::writeVECT(vname = cabins_vect_name, v.in.ogr_flags = "overwrite")

# rasterize with count, creating a new temp_vector
cabins_count_name <- util_v2rast_count_grass(cabins_vect_name, output = "cabins_count", quiet = F, overwrite = T)

# rasterize with count, withour creating a temporary vector
cabins_count_name <- util_v2rast_count_grass(cabins_vect_name, output = "cabins_count",
                                             column = "value", quiet = F, overwrite = T)

# visualize
rgrass7::use_sp()
rgrass7::readRAST(cabins_count_name) %>%
  raster::raster() %>%
  terra::rast() %>%
  plot(main = "Number of private cabins")

# remove rasters created
to_remove_vect <- c(test_region_name, cabins_vect_name)
to_remove_rast <- c(cabins_count_name)
# rgrass7::execGRASS("g.remove", type = "vect", name = to_remove_vect, flags = "f")
# rgrass7::execGRASS("g.remove", type = "rast", name = to_remove_rast, flags = "f")
