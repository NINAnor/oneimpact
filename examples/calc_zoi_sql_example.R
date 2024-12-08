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
DBI::dbExecute(con, "INSTALL spatial from core_nightly; LOAD spatial;")

# write vector of reindeer points to database

# load data in R
data("reindeer")
# register link to data in duckdb
duckdb::duckdb_register(con, "reindeer", reindeer)
# create spatial object in duckdb
DBI::dbExecute(con, "create or replace table reindeer_spat as (select row_number() over () as id, * exclude(x, y), ST_POINT(x,y) as geom from reindeer)")
duckdb::duckdb_unregister(con, "reindeer") # and forget the original dataframe
# check
dplyr::tbl(con, "reindeer_spat")
# add index id and spatial index
DBI::dbExecute(con, "CREATE UNIQUE INDEX reindeer_gid ON reindeer_spat (id);")
DBI::dbExecute(con, "CREATE INDEX reindeer_geometry ON reindeer_spat USING rtree (geom);")

# write vector of cabin points to database - from file
DBI::dbExecute(con, "create or replace table cabins as select * from st_read('inst/vector/reindeer_cabins.gpkg')")
# check
dplyr::tbl(con, "cabins")
# add spatial index
DBI::dbExecute(con, "CREATE INDEX cabins_geom ON cabins USING rtree (geom);")

# compute ZOI of cabins and extract for reindeer points
cum_zoi_cabins <- calc_zoi_sql(con,
                               input_points = "reindeer_spat",
                               infrastructure_layer = "cabins",
                               radius = 5000,
                               type = "bartlett", zoi_metric = "cumulative",
                               input_id = "id",
                               input_geom = "geom", infra_geom = "geom",
                               output_table = NULL,
                               limit = 100,
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
