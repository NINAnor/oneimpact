#' #---
#' Preparing data for the package
#' #---

# Load packages
library(devtools)
library(rgrass7)
library(raster)
library(terra)
library(sf)
library(dplyr)

#------------------
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

# I defined a region_test_influence through the GRASS GUI and using v.in.region output=region_test_influence
rgrass7::execGRASS("g.region", parameters = list(vector = "region_test_influence", align = "private_cabins_rast"),
                   flags = "p")

#------------------
# Study area

# read vector defining study area from GRASS
use_sf()
# study_area <- rgrass7::readVECT(c("region_test_oneimpact_pkg"))
# save externally
sf::st_write(study_area, dsn = "inst/vector/study_area.gpkg", delete_dsn = TRUE)
# re-open
study_area <- sf::st_read("inst/vector/study_area.gpkg")
study_area

# save
# usethis::use_data(study_area, overwrite = TRUE)

# put that data back to GRASS GIS
region_test_oneimpact_pkg <- "region_test_oneimpact_pkg"
rgrass7::writeVECT(study_area, region_test_oneimpact_pkg,
                   v.in.ogr_flags = "overwrite")

#------------------
# Get cabins (already loaded) for this specific study area
cabins_vect_name <- "private_cabins@p_prodchange_envpoints"
rgrass7::execGRASS("v.select", ainput = cabins_vect_name,
                   binput = region_test_oneimpact_pkg, operator = "within",
                   output = "private_cabins_vect", flags = "overwrite")

# get to R
cabins_vect <- rgrass7::readVECT("private_cabins_vect") %>%
  dplyr::select(cat, byggtyp_nbr, kommune, value)
cabins_vect

# save externally
sf::st_write(cabins_vect, dsn = "inst/vector/cabins_vect.gpkg", delete_dsn = TRUE)

#------------------
# Get cabins --- re do here based on the rasterization of the vector above
rgrass7::use_sp()
cabins <- readRAST(c("private_cabins_rast_sub")) %>%
  raster::raster() %>%
  terra::rast()
names(cabins) <- "cabins"
# plot
terra::plot(cabins, col = "black")

# save externally
terra::writeRaster(cabins, filename = "inst/raster/cabins.tif")
# re-open
cabins <- terra::rast("inst/raster/cabins.tif")
cabins

# save
# usethis::use_data(cabins, overwrite = TRUE)

