# prepare reindeer movement data and associated environmental information

# packages
library(readr)
library(lubridate)
library(amt)

# data downloaded from:
# https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.12449
# paper: https://doi.org/10.1111/1365-2656.12449
# datasets: https://doi.org/10.5061/dryad.rg0v3

# read movement data
# coordinates on ETRS 32N crs - EPSG:25832
file_reindeer <- "https://datadryad.org/stash/downloads/file_stream/91560"

reindeer_raw <- readr::read_delim(file_reindeer, delim = "\t",
                                  col_types = c("c", "d", "d", "d", "c", "d", "c", "c")) |>
  dplyr::mutate(sex = "f",
                acquisition_time = lubridate::dmy_hms(acquisition_time)) |>
  dplyr::select(original_animal_id, animal_year_id, sex, utm_x, utm_y, acquisition_time)

# reproject data to UTM 33N to match with enviornmental data
reindeer <- amt::make_track(reindeer_raw, utm_x, utm_y, acquisition_time,
                            crs = sp::CRS("+init=epsg:25832"), all_cols = TRUE) |>
  amt::transform_coords(crs_to = sp::CRS("+init=epsg:25833")) |>
  dplyr::rename(x = x_, y = y_, t = t_)

# save
usethis::use_data(reindeer, overwrite = TRUE)
# now we document by hand

#---------------------
# Get data on the reindeer area
library(sf)
library(terra)
library(DBI)
library(NinaR)

# connect to NINA PostGIS database

source("~/.pgpass")

NinaR::postgreSQLConnect(
  host = "gisdata-db.nina.no",
  dbname = "gisdata",
  username = pg_username,
  password = pg_password
)
rm(pg_username, pg_password)

#--------------
# Get the management area for Austhei - Norway - to match with GPS data
reindeer_area <- sf::st_read(con, DBI::Id(schema = "sam_wrein_ancillary", table = "reindeer_areas"),
                             geometry_column = "geom_e33") |>
  dplyr::select(reindeer_areas_id, name_area, species, geom = geom_e33) |>
  dplyr::filter(name_area == "Setesdal Austhei")

reindeer_area
plot(reindeer_area)

# save
sf::st_write(reindeer_area, dsn = "inst/vector/reindeer_area.shp", delete_dsn = TRUE)
# now we document by hand

# try to read
# cabins count
(s <- system.file("vector/reindeer_area.gpkg", package = "oneimpact"))
v <- terra::vect(s)
terra::plot(v)

#---------------------
# Get maps of environmental covariates

#--------
# get private cabins
cabins_all <- sf::st_read(con, DBI::Id(schema = "sam_env", table = "cabins_private_n50_no"))
# cut
which_within <- which(sf::st_within(cabins_all, reindeer_area, sparse = FALSE))
cabins <- cabins_all |>
  dplyr::slice(which_within) |>
  dplyr::mutate(value = 1) |>
  dplyr::select(gid, buildtype = byggtyp_nbr, city = kommune, value)

cabins
plot(cabins)

# save
sf::st_write(cabins, dsn = "inst/vector/reindeer_cabins.gpkg", delete_dsn = TRUE)
# cabins <- sf::st_read(dsn = "inst/vector/reindeer_cabins.gpkg")
# now we document by hand

# try to read
# cabins count
(s <- system.file("vector/reindeer_cabins.gpkg", package = "oneimpact"))
v <- terra::vect(s)
terra::plot(v)

#--------
# get public roads
roads_public_all <- sf::st_read(con, DBI::Id(schema = "sam_env", table = "roads_public_renrein_no"))
# cut
roads_public <- sf::st_intersection(roads_public_all, reindeer_area) |>
  dplyr::mutate(value = 1) |>
  dplyr::select(id, name = gatenavn, publ_priv, traffic_bin, name_area, value)

roads_public
plot(roads_public[4])

# save
sf::st_write(roads_public, dsn = "inst/vector/reindeer_roads_public.gpkg", delete_dsn = TRUE)
# roads_public <- sf::st_read(dsn = "inst/vector/reindeer_roads_public.gpkg")
# now we document by hand

# try to read
# cabins count
(s <- system.file("vector/reindeer_roads_public.gpkg", package = "oneimpact"))
v <- sf::st_read(s)
plot(v[4])

#--------
# get private roads
roads_private_all <- sf::st_read(con, DBI::Id(schema = "sam_env", table = "roads_private_renrein_no"))
# cut
roads_private <- sf::st_intersection(roads_private_all, reindeer_area) |>
  dplyr::mutate(value = 1) |>
  dplyr::select(id, name = gatenavn, publ_priv, traffic_bin, name_area, value)

roads_private
plot(roads_private[4])

# save
sf::st_write(roads_private, dsn = "inst/vector/reindeer_roads_private.gpkg", delete_dsn = TRUE)
roads_private <- sf::st_write(dsn = "inst/vector/reindeer_roads_private.gpkg")
# now we document by hand

# try to read
# cabins count
(s <- system.file("vector/roads_private.gpkg", package = "oneimpact"))
v <- sf::st_read(s)
plot(v[4])

#---------------------------------
# raster data

# base raster
(s <- system.file("vector/reindeer_area.gpkg", package = "oneimpact"))
v <- terra::vect(s)
terra::plot(v)

base_rast <- terra::rast(v, res = 100)

#---
# raster - cabins presence
(s <- system.file("vector/cabins.gpkg", package = "oneimpact"))
v <- terra::vect(s)

# cabins_presence <- terra::rasterize(v, base_rast, fun = length)
# cabins_presence[is.na(cabins_presence)] <- 0
# cabins_presence[cabins_presence > 0] <- 1
cabins_presence <- terra::rasterize(v, base_rast, field = "value")
plot(cabins_presence)

cabins_count <- terra::rasterize(v, base_rast, fun = length)
cabins_count[is.na(cabins_count)] <- 0
plot(cabins_count)

# save externally
terra::writeRaster(cabins_count, filename = "inst/raster/cabins_count.tif")
# re-open
cabins_count <- terra::rast("inst/raster/cabins_count.tif")
cabins_count

library(oneimpact)
calc_zoi_nearest(cabins_presence, type = "euclidean")
