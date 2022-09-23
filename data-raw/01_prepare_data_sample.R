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
rgrass7::execGRASS("g.region", vector = "region_test_oneimpact_pkg", res = "100",
                   flags = "p")

#------------------
# Study area

# read vector defining study area from GRASS
use_sf()
sample_area <- rgrass7::read_VECT(c("region_test_oneimpact_pkg"))
# save externally
sf::st_write(sample_area, dsn = "inst/vector/sample_area.gpkg", delete_dsn = TRUE)
# re-open
sample_area <- sf::st_read("inst/vector/sample_area.gpkg")
sample_area

# save
# usethis::use_data(study_area, overwrite = TRUE)

# put that data back to GRASS GIS
# region_test_oneimpact_pkg <- "region_test_oneimpact_pkg"
# rgrass7::writeVECT(sample_area, region_test_oneimpact_pkg,
#                    v.in.ogr_flags = "overwrite")

#------------------
# Get cabins (already loaded) for this specific study area
cabins_vect_name <- "private_cabins@p_prodchange_envpoints"
rgrass7::execGRASS("v.select", ainput = cabins_vect_name,
                   binput = region_test_oneimpact_pkg, operator = "within",
                   output = "private_cabins_sample_area", flags = "overwrite")

# get to R
cabins_vect <- rgrass7::read_VECT("private_cabins_sample_area") |>
  dplyr::select(cat, buildtype = byggtyp_nbr, city = kommune, value)
cabins_vect

# save externally
sf::st_write(cabins_vect, dsn = "inst/vector/sample_area_cabins.gpkg", delete_dsn = TRUE)

#------------------
# Get cabins

#------
# cabins 1/null

# load vector
cabins_vect <- terra::vect("inst/vector/sample_area_cabins.gpkg")
cabins_vect
# load old raster
r <- terra::rast("inst/raster/cabins_sample.tif")
# rasterize
cabins_rast <- terra::rasterize(cabins_vect, r, field = "value")
# plot
terra::plot(cabins_rast, col = "black")

# save externally
terra::writeRaster(cabins_rast, filename = "inst/raster/sample_area_cabins.tif")
# re-open
cabins <- terra::rast("inst/raster/sample_area_cabins.tif")
cabins

# save
# usethis::use_data(cabins, overwrite = TRUE)

#------
# cabins count

# load vector
cabins_vect <- terra::vect("inst/vector/sample_area_cabins.gpkg")
cabins_vect
# load raster
r <- terra::rast("inst/raster/sample_area_cabins.tif")

cabins_count <- terra::rasterize(cabins_vect, r, fun = length)
cabins_count[is.na(cabins_count)] <- 0
plot(cabins_count)

# save externally
terra::writeRaster(cabins_count, filename = "inst/raster/sample_area_cabins_count.tif")
# re-open
cabins_count <- terra::rast("inst/raster/sample_area_cabins_count.tif")
cabins_count

#------------------
# Get roads public

#--------
# Vector - not included in the package, just for checking here
roads_vect_name <- "roads_public_no@p_sam_transport_urban"
rgrass7::execGRASS("v.select", ainput = roads_vect_name,
                   binput = "region_test_oneimpact_pkg", operator = "within",
                   output = "roads_public_sample_area", flags = "overwrite")

# get to R
roads_public_vect <- rgrass7::read_VECT("roads_public_sample_area") |>
  sf::st_as_sf() |>
  dplyr::mutate(value = 1) |>
  dplyr::select(id, name = gatenavn, publ_priv, traffic_bin, value)
plot(roads_public_vect[1])

# save externally
sf::st_write(roads_public_vect, dsn = "inst/vector/sample_area_roads.gpkg", delete_dsn = TRUE)

#------
# roads 1/null

# load vector
roads_public_vect <- terra::vect("inst/vector/sample_area_roads.gpkg")
roads_public_vect
# load cabins raster
r <- terra::rast("inst/raster/sample_area_cabins.tif")
# rasterize
roads_rast <- terra::rasterize(roads_public_vect, r, field = "value")
# plot
terra::plot(roads_rast, col = "black")

# save externally
terra::writeRaster(roads_rast, filename = "inst/raster/sample_area_roads.tif")
# re-open
cabins <- terra::rast("inst/raster/sample_area_roads.tif")
cabins

# save
# usethis::use_data(cabins, overwrite = TRUE)

#------
# cabins count

# load vector
cabins_vect <- terra::vect("inst/vector/sample_area_cabins.gpkg")
cabins_vect
# load raster
r <- terra::rast("inst/raster/sample_area_cabins.tif")

cabins_count <- terra::rasterize(cabins_vect, r, fun = length)
cabins_count[is.na(cabins_count)] <- 0
plot(cabins_count)

# save externally
terra::writeRaster(cabins_count, filename = "inst/raster/sample_area_cabins_count.tif")
# re-open
cabins_count <- terra::rast("inst/raster/sample_area_cabins_count.tif")
cabins_count


#------------------
# Get roads private

#--------
# Vector - # not included in the package, but it could be included

roads_priv_vect_name <- "roads_private_no@p_sam_transport_urban"
rgrass7::execGRASS("v.select", ainput = roads_priv_vect_name,
                   binput = "region_test_oneimpact_pkg", operator = "within",
                   output = "roads_private_sample_area", flags = "overwrite")

# get to R
roads_private_vect <- rgrass7::read_VECT("roads_private_sample_area") |>
  sf::st_as_sf() |>
  dplyr::mutate(value = 1) |>
  dplyr::select(id, name = gatenavn, publ_priv, traffic_bin, value)
plot(roads_private_vect[1])
# not included in the package, but it could be included


#------------------------
# test
library(oneimpact)

(f <- system.file("raster/sample_area_roads.tif", package = "oneimpact"))
roads <- terra::rast(f)

rn <- calc_zoi_nearest(roads, type = "euclidean")
plot(rn)

rc <- calc_zoi_cumulative(roads, type = "bartlett", radius = 1000,
                          zeroAsNA = TRUE)
plot(rc)
