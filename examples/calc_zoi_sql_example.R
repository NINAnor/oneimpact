library(oneimpact)

library(DBI)
library(dplyr)
library(sf)
library(terra)
# install.packages("duckdb")
library(duckdb)

#---
# set up connection and files

# connection - in memory
con <- DBI::dbConnect(duckdb())
DBI::dbExecute(con, "INSTALL spatial; LOAD spatial;")
# source("~/.pgpass")
# NinaR::postgreSQLConnect(
#   host = "gisdata-db.nina.no",
#   dbname = "gisdata",
#   username = pg_username,
#   password = pg_password
# )
# rm(pg_username, pg_password)

# write vector of reindeer points to database
data("reindeer")
rein_spat <- sf::st_as_sf(reindeer, coords = c("x", "y"), crs = 25833) |>
  dplyr::mutate(gid = 1:nrow(reindeer))
DBI::dbWriteTable(con, "reindeer", rein_spat, overwrite = TRUE)
# check
dplyr::tbl(con, "reindeer")
# add index gid and spatial index
DBI::dbExecute(con, "CREATE UNIQUE INDEX reindeer_gid ON reindeer (gid);")
DBI::dbExecute(con, "CREATE INDEX reindeer_geometry ON reindeer USING GIST (geometry);")

# write vector of cabin points to database
(f <- system.file("vector/reindeer_cabins.gpkg", package = "oneimpact"))
cabins <- sf::st_read(f)
DBI::dbWriteTable(con, "cabins", cabins)
# check
dplyr::tbl(con, "cabins")
# add spatial index
DBI::dbExecute(con, "CREATE INDEX cabins_geom ON cabins USING GIST (geom);")

# compute ZOI of cabins and extract for reindeer points
cum_zoi_cabins <- calc_zoi_sql(input_points = "reindeer",
                               infrastructure_layer = "cabins",
                               radius = 5000,
                               type = "bartlett", zoi_metric = "cumulative",
                               input_id = "gid",
                               input_geom = "geometry", infra_geom = "geom",
                               output_table = NULL,
                               #limit = 100,
                               verbose = TRUE)
cum_zoi_cabins

#-----------------
# compare to the raster approach

# compute ZOI
f <- system.file("vector/reindeer_cabins.gpkg", package = "oneimpact")
cabins <- terra::vect(f)
rr <- terra::rast(xmin = terra::ext(cabins)[1], resolution = 100,
                  extent = terra::ext(cabins), crs = terra::crs(cabins))
cabins_count <- terra::rasterize(cabins, rr, fun = length)
cabins_count <- terra::ifel(is.na(cabins_count), 0, cabins_count)
plot(cabins_count)

cumzoi_linear <- calc_zoi_cumulative(cabins_count, type = "bartlett", radius = 5000)
plot(cumzoi_linear)

# extract
reindeer_cabins <- terra::extract(cumzoi_linear, terra::vect(rein_spat))

plot(cumzoi_linear)
plot(terra::vect(rein_spat), add = T)

# approximately the same
cbind(reindeer_cabins, dplyr::arrange(cum_zoi_cabins, gid))
