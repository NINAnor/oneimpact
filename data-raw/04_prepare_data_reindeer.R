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
                            crs = 25832, all_cols = TRUE) |>
  amt::transform_coords(crs_to = 25833) |>
  dplyr::rename(x = x_, y = y_, t = t_) |>
  dplyr::arrange(original_animal_id, t)

reindeer <- tibble::as_tibble(reindeer)
class(reindeer)

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


#------------------------------------
# Annotated reindeer data for SSF
library(terra)
library(amt)
library(dplyr)
library(tidyr)
data(reindeer)

reindeer |>
  dplyr::filter(animal_year_id == 33580) |>
  amt::make_track(x, y, t, crs = 25833, all_cols = TRUE) |>
  # dplyr::rename(x_ = x, y_ = y, t_ = t) |>
  dplyr::filter(lubridate::month(t_) == 7) |>
  amt::track_resample(rate = hours(3), tolerance = hours(1)) |>
  amt::filter_min_n_burst() |>
  amt::steps_by_burst(keep_cols = "start") |>
  amt::random_steps()

# use-availability setup for ssf
r <- reindeer |>
  amt::make_track(x, y, t, crs = 25833, all_cols = TRUE) |>
  # dplyr::rename(x_ = x, y_ = y, t_ = t) |>
  dplyr::filter(lubridate::month(t_) == 7) |> # only July
  tidyr::nest(.by = animal_year_id) |>
  dplyr::mutate(ua = purrr::map(data, function(x)
    x |>
      amt::track_resample(rate = hours(3), tolerance = hours(1)) |>
      amt::filter_min_n_burst() |>
      amt::steps_by_burst(keep_cols = "start") |>
      amt::random_steps())) |>
  dplyr::select(-data) |>
  tidyr::unnest(ua) |>
  dplyr::mutate(step_id = paste(animal_year_id, burst_, step_id_, sep = "_"))

# compute zoi
reindeer_area_v <- terra::vect(system.file("vector/reindeer_area.gpkg", package = "oneimpact"))
# extent <- terra::buffer(reindeer_area_v, max(r$sl_)*1.1) |>
#   terra::ext()
extent <- r |>
  terra::vect(geom = c("x1_", "y1_")) |>
  terra::buffer(max(r$sl_)*1.1) |>
  terra::ext()
ress <- 200
nrows <- floor(range(extent)[2]/ress+1)
ncols <- floor(range(extent)[1]/ress+1)
rr <- terra::rast(extent, nrows = nrows, ncols = ncols)
# cabins
reindeer_cabins_v <- terra::vect(system.file("vector/reindeer_cabins.gpkg", package = "oneimpact"))
reindeer_cabins <- terra::rasterize(reindeer_cabins_v, rr, fun = length)
cabins_zoi <- calc_zoi_cumulative(reindeer_cabins, radius = c(seq(1000, 3000, 1000)),
                                  type = "Gauss", zeroAsNA = TRUE)
names(cabins_zoi) <- paste0("cabins", seq(1000, 3000, 1000))
terra::plot(cabins_zoi)
# roads
reindeer_roads_v <- terra::vect(system.file("vector/reindeer_roads_public.gpkg", package = "oneimpact"))
reindeer_roads <- terra::rasterize(reindeer_roads_v, rr, field = "value")
roads_zoi <- calc_zoi_cumulative(reindeer_roads, radius = c(seq(1000, 3000, 1000)),
                                 type = "Gauss", zeroAsNA = TRUE)
names(roads_zoi) <- paste0("roads", seq(1000, 3000, 1000))
plot(roads_zoi)

# extract
st <- terra::extract(c(cabins_zoi, roads_zoi), terra::vect(r, geom = c("x1_", "y1_")))
names(st)[-1] <- paste0("start_", names(st)[-1])
end <- terra::extract(c(cabins_zoi, roads_zoi), terra::vect(r, geom = c("x2_", "y2_")))
names(end)[-1] <- paste0("end_", names(end)[-1])

# anottated data
data_annotated <- dplyr::bind_cols(r, st[-1], end[-1])
data_annotated[grep("cabins|roads", names(data_annotated))[1]:ncol(data_annotated)] <-
  apply(data_annotated[grep("cabins|roads", names(data_annotated))[1]:ncol(data_annotated)], 2, scale)
data_annotated$case_ <- ifelse(data_annotated$case_, 1, 0)

# plot
data_annotated_v <- terra::vect(data_annotated, geom = c("x2_", "y2_"))
plot(data_annotated_v[data_annotated_v$case_ == F], col = "grey", cex = 0.3)
points(data_annotated_v[data_annotated_v$case_ == T], col = "red", cex = 0.3)

# save
reindeer_ssf <- data_annotated
usethis::use_data(reindeer_ssf, overwrite = TRUE)
# now we document by hand

#------------------------------------
# Annotated reindeer data for SSF
library(terra)
library(amt)
library(dplyr)
library(tidyr)
data(reindeer)

# use-availability setup for ssf
r <- reindeer |>
  amt::make_track(x, y, t, crs = 25833, all_cols = TRUE) |>
  # dplyr::rename(x_ = x, y_ = y, t_ = t) |>
  dplyr::filter(lubridate::month(t_) == 7) |> # only July
  tidyr::nest(.by = animal_year_id) |>
  dplyr::mutate(ua = purrr::map(data, function(x)
    x |>
      amt::track_resample(rate = hours(3), tolerance = hours(1)) |>
      amt::filter_min_n_burst() |>
      amt::steps_by_burst(keep_cols = "start") |>
      amt::random_steps())) |>
  dplyr::select(-data) |>
  tidyr::unnest(ua) |>
  dplyr::mutate(step_id = paste(animal_year_id, burst_, step_id_, sep = "_"))

# compute zoi
reindeer_area_v <- terra::vect(system.file("vector/reindeer_area.gpkg", package = "oneimpact"))
# extent <- terra::buffer(reindeer_area_v, max(r$sl_)*1.1) |>
#   terra::ext()
extent <- r |>
  terra::vect(geom = c("x1_", "y1_")) |>
  terra::buffer(max(r$sl_)*1.1) |>
  terra::ext()
ress <- 200
nrows <- floor(range(extent)[2]/ress+1)
ncols <- floor(range(extent)[1]/ress+1)
rr <- terra::rast(extent, nrows = nrows, ncols = ncols)
# cabins
reindeer_cabins_v <- terra::vect(system.file("vector/reindeer_cabins.gpkg", package = "oneimpact"))
reindeer_cabins <- terra::rasterize(reindeer_cabins_v, rr, fun = length)
cabins_zoi <- calc_zoi_cumulative(reindeer_cabins, radius = c(seq(1000, 3000, 1000)),
                                  type = "Gauss", zeroAsNA = TRUE)
names(cabins_zoi) <- paste0("cabins", seq(1000, 3000, 1000))
terra::plot(cabins_zoi)
# roads
reindeer_roads_v <- terra::vect(system.file("vector/reindeer_roads_public.gpkg", package = "oneimpact"))
reindeer_roads <- terra::rasterize(reindeer_roads_v, rr, field = "value")
roads_zoi <- calc_zoi_cumulative(reindeer_roads, radius = c(seq(1000, 3000, 1000)),
                                 type = "Gauss", zeroAsNA = TRUE)
names(roads_zoi) <- paste0("roads", seq(1000, 3000, 1000))
plot(roads_zoi)

# extract
st <- terra::extract(c(cabins_zoi, roads_zoi), terra::vect(r, geom = c("x1_", "y1_")))
names(st)[-1] <- paste0("start_", names(st)[-1])
end <- terra::extract(c(cabins_zoi, roads_zoi), terra::vect(r, geom = c("x2_", "y2_")))
names(end)[-1] <- paste0("end_", names(end)[-1])

# anottated data
data_annotated <- dplyr::bind_cols(r, st[-1], end[-1])
data_annotated[grep("cabins|roads", names(data_annotated))[1]:ncol(data_annotated)] <-
  apply(data_annotated[grep("cabins|roads", names(data_annotated))[1]:ncol(data_annotated)], 2, scale)
data_annotated$case_ <- ifelse(data_annotated$case_, 1, 0)

# plot
data_annotated_v <- terra::vect(data_annotated, geom = c("x2_", "y2_"))
plot(data_annotated_v[data_annotated_v$case_ == F], col = "grey", cex = 0.3)
points(data_annotated_v[data_annotated_v$case_ == T], col = "red", cex = 0.3)

# save
reindeer_ssf <- data_annotated
usethis::use_data(reindeer_ssf, overwrite = TRUE)
# now we document by hand

#------------------------------------
# Annotated reindeer data for RSF
# packages
library(oneimpact)
library(glmnet)
library(terra)
library(amt)
data(reindeer)

# use-availability setup for rsf
r <- reindeer |>
  amt::make_track(x, y, t, crs = 25833, all_cols = TRUE) |>
  # dplyr::rename(x_ = x, y_ = y, t_ = t) |>
  dplyr::filter(lubridate::month(t_) == 7)

hr_buff1km <- reindeer |>
  amt::make_track(x, y, t, crs = 25833, all_cols = TRUE) |>
  amt::hr_mcp(level = 1)
hr_buff1km$mcp <- sf::st_buffer(hr_buff1km$mcp, 1000)
rp <- amt::random_points(hr_buff1km, n = nrow(r) * 10)

# dataset
r <- r |>
  dplyr::mutate(case_ = TRUE) |>
  dplyr::select(case_, x_, y_) |>
  dplyr::bind_rows(rp)

# compute zoi
# reindeer_area_v <- terra::vect(system.file("vector/reindeer_area.gpkg", package = "oneimpact"))
# extent <- terra::buffer(reindeer_area_v, max(r$sl_)*1.1) |>
#   terra::ext()
extent <- r |>
  terra::vect(geom = c("x_", "y_")) |>
  terra::buffer(3000) |>
  terra::ext()
ress <- 200
nrows <- floor(range(extent)[2]/ress+1)
ncols <- floor(range(extent)[1]/ress+1)
rr <- terra::rast(extent, nrows = nrows, ncols = ncols)
# cabins
reindeer_cabins_v <- terra::vect(system.file("vector/reindeer_cabins.gpkg", package = "oneimpact"))
reindeer_cabins <- terra::rasterize(reindeer_cabins_v, rr, fun = length)
cabins_zoi <- calc_zoi_cumulative(reindeer_cabins, radius = c(seq(1000, 3000, 1000)),
                                  type = "Gauss", zeroAsNA = TRUE)
names(cabins_zoi) <- paste0("cabins", seq(1000, 3000, 1000))
terra::plot(cabins_zoi)
# roads
reindeer_roads_v <- terra::vect(system.file("vector/reindeer_roads_public.gpkg", package = "oneimpact"))
reindeer_roads <- terra::rasterize(reindeer_roads_v, rr, field = "value")
roads_zoi <- calc_zoi_cumulative(reindeer_roads, radius = c(seq(1000, 3000, 1000)),
                                 type = "Gauss", zeroAsNA = TRUE)
names(roads_zoi) <- paste0("roads", seq(1000, 3000, 1000))
plot(roads_zoi)

# extract
pts <- terra::extract(c(cabins_zoi, roads_zoi), terra::vect(r, geom = c("x_", "y_")))

# anottated data
reindeer_rsf <- dplyr::bind_cols(r, pts[-1])
reindeer_rsf[grep("cabins|roads", names(reindeer_rsf))[1]:ncol(reindeer_rsf)] <-
  apply(reindeer_rsf[grep("cabins|roads", names(reindeer_rsf))[1]:ncol(reindeer_rsf)], 2, scale)
reindeer_rsf$case_ <- ifelse(reindeer_rsf$case_, 1, 0)

# plot
reindeer_rsf_v <- terra::vect(reindeer_rsf, geom = c("x_", "y_"))
plot(reindeer_rsf_v[reindeer_rsf$case_ == F], col = "grey", cex = 0.3)
points(reindeer_rsf_v[reindeer_rsf$case_ == T], col = "red", cex = 0.3)

# compute squared values for some columns
reindeer_rsf$norway_pca_klima_axis1_sq <- reindeer_rsf$norway_pca_klima_axis1**2
reindeer_rsf$norway_pca_klima_axis2_sq <- reindeer_rsf$norway_pca_klima_axis2**2

# save
usethis::use_data(reindeer_rsf, overwrite = TRUE)
# now we document by hand

#---------------
# adding data for RSF in Hardangervidda from Niebuhr et al 2023 MEE

# One pop, pts
load("/data/P-Prosjekter/41203800_oneimpact/04_tools/support_oneimpact/cuminf_zoi_GPS_dataset_annotated.rda")
dat
str(dat)

#------------------------------
# One fit - best fit from Niebuhr et al 2023

# table of fits
load(file = "/data/P-Prosjekter/41203800_oneimpact/04_tools/support_oneimpact/cuminf_zoi_results_rsf_priv_pub_cabins_glm.rda")
head(multi_infra_model_comparison_df)
multi_infra_best_model_call
f <- as.formula(multi_infra_best_model_call)

# formula with all variables
f <- as.formula(multi_infra_best_model_call)
f <- as.character(f) |> gsub(pattern = "cumulative_threshold_10000|cumulative_exp_decay_20000", replacement = "XXX")
f <- as.formula(paste(f[2], f[1], f[3]))
zois <- c(100, 250, 500, 1000, 2500, 5000, 10000, 20000)
ff <- add_zoi_formula(f, zoi_radius = zois, pattern = "XXX", type = c("cumulative_exp_decay", "nearest_exp_decay"),
                      separator = "_", grid = TRUE)

f <- ff$formula
grid_zoi <- ff$grid

# get data
reindeer_rsf <- dat[,all.vars(f)]
# save
usethis::use_data(reindeer_rsf, overwrite = TRUE)
# now we document by hand

#------------------------------
# all ZOI covariates from Niebuhr et al 2023


# table of fits
load(file = "/data/P-Prosjekter/41203800_oneimpact/04_tools/support_oneimpact/cuminf_zoi_results_rsf_priv_pub_cabins_glm.rda")
head(multi_infra_model_comparison_df)
multi_infra_best_model_call
f <- as.formula(multi_infra_best_model_call)


# formula with all variables
f <- as.formula(multi_infra_best_model_call)
f <- as.character(f) |> gsub(pattern = "cumulative_threshold_10000|cumulative_exp_decay_20000", replacement = "XXX")
f <- as.formula(paste(f[2], f[1], f[3]))
zois <- c(100, 250, 500, 1000, 2500, 5000, 10000, 20000)
ff <- add_zoi_formula(f, zoi_radius = zois, pattern = "XXX", type = c("cumulative_exp_decay", "nearest_exp_decay"),
                      separator = "_", grid = TRUE)

f <- ff$formula
grid_zoi <- ff$grid

# get data
reindeer_rsf <- dat[,all.vars(f)]
# save

#-----------------------------------------------------

# add linear infrastructure to Hardanger RSF example
library(terra)

# prepared data
reindeer_rsf
# saveRDS(reindeer_rsf, "data-raw/reindeer_rsf_old_only_cabins_pca_landuse_hardanger.rda")

# base data
load("/data/P-Prosjekter/41203800_oneimpact/04_tools/support_oneimpact/cuminf_zoi_GPS_dataset_annotated.rda")
dat
str(dat)

# check
nrow(reindeer_rsf)
nrow(dat)

# subset
dat <- dat |>
  dplyr::select(points_id, x33, y33, herd, use)

# annotate
r1 <- terra::rast("data-raw/rast_predictors_hardanger_100.tif")
r2 <- terra::rast("data-raw/rast_predictors_hardanger_100_linear.tif")
names(r2) <- names(r2) |>
  sub(pattern = "bin_", replacement = "") |>
  sub(pattern = "inf_", replacement = "")
ext1 <- terra::extract(r1, dat[, c("x33", "y33")])
ext2 <- terra::extract(r2, dat[, c("x33", "y33")])

colnames(reindeer_rsf) %in% colnames(ext1)
colnames(ext1) %in% colnames(reindeer_rsf)

# combine
names(ext1)
ext1[, c(2:5, 39:40, 6, 6+c(1,5,7,2,6,8,3,4), 14+c(1,5,7,2,6,8,3,4), 22+c(1,5,7,2,6,8,3,4), 30+c(1,5,7,2,6,8,3,4))] |> names()
ext1 <- ext1[, c(2:5, 39:40, 6, 6+c(1,5,7,2,6,8,3,4), 14+c(1,5,7,2,6,8,3,4), 22+c(1,5,7,2,6,8,3,4), 30+c(1,5,7,2,6,8,3,4))]

names(ext2)
ext2[, c(1+c(1,5,7,2,6,8,3,4), 9+c(1,5,7,2,6,8,3,4), 17+c(1,5,7,2,6,8,3,4), 25+c(1,5,7,2,6,8,3,4), 33+c(1,5,7,2,6,8,3,4), 41+c(1,5,7,2,6,8,3,4))] |> names()
ext2 <- ext2[, c(1+c(1,5,7,2,6,8,3,4), 9+c(1,5,7,2,6,8,3,4), 17+c(1,5,7,2,6,8,3,4), 25+c(1,5,7,2,6,8,3,4), 33+c(1,5,7,2,6,8,3,4), 41+c(1,5,7,2,6,8,3,4))]

reindeer_rsf <- cbind(dat[, 5, drop = FALSE], ext1, ext2)
names(reindeer_rsf)

names(reindeer_rsf) <- names(reindeer_rsf) |>
  sub(pattern = "decay_", replacement = "decay")

# export for the package
usethis::use_data(reindeer_rsf, overwrite = TRUE)


##-----------------------
# for Ron

library(arrow)
dat_ron <- dat |>
  dplyr::select(points_id, x33, y33, herd, use)
dat_ron |>
  arrow::write_parquet("/data/P-Prosjekter/41203800_oneimpact/02_sam/02_data/bio/processed/data_hardanger_ron.parquet")

rr <- terra::rast("data-raw/rast_predictors_hardanger_100.tif")
rr[[1:5]] |>
  terra::writeRaster("/data/P-Prosjekter/41203800_oneimpact/02_sam/02_data/bio/processed/env_pcas_landcover.tif")

##-----------------------
# Adding vectors for Hardanger with a buffer of 20 km around the area

library(sf)
library(terra)
library(DBI)
library(NinaR)
library(rgrass)

# connect to NINA PostGIS database

source("~/.pgpass")

NinaR::postgreSQLConnect(
  host = "gisdata-db.nina.no",
  dbname = "gisdata",
  username = pg_username,
  password = pg_password
)
rm(pg_username, pg_password)

ms <- "u_bb_cuminf"
NinaR::grassConnect(mapset = ms)

#--------------
# Set region

rgrass::execGRASS("r.mask", flags = "r")

# region
rgrass::execGRASS("v.info", map = "reindeer_areas_no_2023@p_sam_tools")

rgrass::execGRASS("v.extract",
                  input = "reindeer_areas_no_2023@p_sam_tools",
                  output = "hardanger_temp_vector",
                  where = "name_area=\"'Hardangervidda'\"",
                  flags = "overwrite")

(f <- system.file("raster/hardanger_rast_predictors_500m.tif", package = "oneimpact"))
r <- terra::rast(f)
rgrass::write_RAST(r[[1]], "hardanger_raster_ref", flags = c("o", "overwrite"))

rgrass::execGRASS("g.region", raster = "raster_ref.1", flags = c("a", "p"), res = "100")
rgrass::execGRASS("g.region", n = "n+12000", e = "e+12000", s = "s-12000", w = "w-12000", flags = c("a", "p"))

# get dimensions

#--------------
# Get the management area for Austhei - Norway - to match with GPS data
reindeer_area_original <- sf::st_read(con, DBI::Id(schema = "sam_wrein_ancillary", table = "reindeer_areas"),
                             geometry_column = "geom_e33") |>
  dplyr::select(reindeer_areas_id, name_area, species, geom = geom_e33) |>
  dplyr::filter(name_area == "Hardangervidda")

reindeer_area <- reindeer_area_original |>
  sf::st_buffer(dist = 12000)

reindeer_area
plot(reindeer_area[1])
plot(reindeer_area_original, add = T)

sf::st_write(reindeer_area_original[,-3], dsn = "inst/vector/hardanger_polygon.gpkg", delete_dsn = TRUE)
# cabins <- sf::st_read(dsn = "inst/vector/hardanger_cabins_private.gpkg")

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
plot(terra::vect(cabins[1]))

# save
sf::st_write(cabins, dsn = "inst/vector/hardanger_cabins_private.gpkg", delete_dsn = TRUE)
# cabins <- sf::st_read(dsn = "inst/vector/hardanger_cabins_private.gpkg")
# now we document by hand

# try to read
# cabins count
(s <- system.file("vector/hardanger_cabins_private.gpkg", package = "oneimpact"))
v <- terra::vect(s)
terra::plot(v)

# compute ZOI
cabins_grass <- "hardanger_buff_12km_cabins_private"
rgrass::write_VECT(terra::vect(cabins), cabins_grass)
oneimpact::grass_v2rast_count(x = cabins_grass, overwrite = TRUE)

xx <- paste0(cabins_grass, "_count")
rad <- c(100, 250, 500, 1000, 2500, 5000, 10000)
cum_layers <- calc_zoi_cumulative(x = xx,
                    radius = rad,
                    type = "exp_decay",
                    where = "GRASS",
                    g_overwrite = TRUE,
                    verbose = TRUE)
cum_layers <- paste0("hardanger_buff_12km_cabins_private_count_zoi_cumulative_exp_decay", rad)

xx2 <- sub(pattern = "_count", replacement = "_1null", xx)
rgrass::execGRASS("r.mapcalc", expression = paste0(xx2, " = if(", xx, " == 0, null(), 1)"), flags = "overwrite")
near_layers <- calc_zoi_nearest(x = xx2,
                               radius = rad,
                               type = "exp_decay",
                               where = "GRASS",
                               g_overwrite = TRUE,
                               verbose = TRUE)
near_layers <- paste0("hardanger_buff_12km_cabins_private_1null_zoi_nearest_exp_decay", rad)

rasters_temp <- rgrass::read_RAST(c(cum_layers, near_layers), return_format = "terra")
plot(rasters_temp[[c(2,6,9,13)]])

#--------
# get public cabins
cabins_all <- sf::st_read(con, DBI::Id(schema = "sam_env", table = "cabins_public_large_no"))
hotels_all <- sf::st_read(con, DBI::Id(schema = "sam_env", table = "hotels_no"))
# cut
which_within <- which(sf::st_within(cabins_all, reindeer_area, sparse = FALSE))
which_within_hotels <- which(sf::st_within(hotels_all, reindeer_area, sparse = FALSE))
cabins <- cabins_all |>
  dplyr::slice(which_within) |>
  dplyr::mutate(value = 1) |>
  dplyr::select(gid, buildtype = byggtyp_nbr, city = kommune, value) |>
  dplyr::bind_rows(
    hotels_all |>
      dplyr::slice(which_within_hotels) |>
      dplyr::mutate(value = 1) |>
      dplyr::select(gid, buildtype = byggtyp_nbr, city = kommune, value)
  )

cabins
plot(terra::vect(cabins))

# save
sf::st_write(cabins, dsn = "inst/vector/hardanger_cabins_public.gpkg", delete_dsn = TRUE)
# cabins <- sf::st_read(dsn = "inst/vector/hardanger_cabins_public.gpkg")
# now we document by hand

# try to read
# cabins count
(s <- system.file("vector/hardanger_cabins_public.gpkg", package = "oneimpact"))
v <- terra::vect(s)
terra::plot(v)

# compute ZOI
cabins_public_grass <- "hardanger_buff_12km_cabins_public_hotels"
rgrass::write_VECT(terra::vect(cabins), cabins_public_grass)
oneimpact::grass_v2rast_count(x = cabins_public_grass, overwrite = TRUE)

xx <- paste0(cabins_public_grass, "_count")
rad <- c(100, 250, 500, 1000, 2500, 5000, 10000)
cum_layers <- calc_zoi_cumulative(x = xx,
                                  radius = rad,
                                  type = "exp_decay",
                                  where = "GRASS",
                                  g_overwrite = TRUE,
                                  verbose = TRUE)
cum_layers <- paste0("hardanger_buff_12km_cabins_public_hotels_count_zoi_cumulative_exp_decay", rad)

xx2 <- sub(pattern = "_count", replacement = "_1null", xx)
rgrass::execGRASS("r.mapcalc", expression = paste0(xx2, " = if(", xx, " == 0, null(), 1)"), flags = "overwrite")
near_layers <- calc_zoi_nearest(x = xx2,
                                radius = rad,
                                type = "exp_decay",
                                where = "GRASS",
                                g_overwrite = TRUE,
                                verbose = TRUE)
near_layers <- paste0("hardanger_buff_12km_cabins_public_hotels_1null_zoi_cumulative_exp_decay", rad)

rasters_temp <- rgrass::read_RAST(c(cum_layers, near_layers), return_format = "terra")
plot(rasters_temp[[c(2,6,9,13)]])

#--------
# get trails
trails_all <- sf::st_read(con, DBI::Id(schema = "sam_env", table = "trails_summer_renrein_no"))
# cut
trails <- sf::st_intersection(trails_all, reindeer_area) |>
  dplyr::mutate(value = 1) |>
  dplyr::select(gid, area, traffic_bin, value)

trails
plot(trails[4])

# save
sf::st_write(trails, dsn = "inst/vector/hardanger_trails.gpkg", delete_dsn = TRUE)
# trails <- sf::st_read(dsn = "inst/vector/hardanger_trails.gpkg")
# now we document by hand

# try to read
# cabins count
(s <- system.file("vector/hardanger_trails.gpkg", package = "oneimpact"))
v <- sf::st_read(s)
plot(v[4])

# compute ZOI
trails_grass <- "hardanger_buff_12km_trails"
rgrass::write_VECT(terra::vect(trails), trails_grass)
rgrass::execGRASS("v.to.rast", input = trails_grass,
                  output = paste0(trails_grass, "_rast_1null"),
                  use = "value", value = 1,
                  flags = c("d", "overwrite"))
rgrass::execGRASS("r.mapcalc",
                  expression = paste0(
                    trails_grass, "_rast_bin = if(isnull(",
                    trails_grass, "_rast_1null), 0, ", trails_grass, "_rast_1null)"),
                  flags = "overwrite")

# rgrass::read_RAST(paste0(trails_grass, "_rast_bin")) |> plot()
# rgrass::read_RAST(paste0(trails_grass, "_rast_1null")) |> plot()

xx <- paste0(trails_grass, "_rast_bin")
rad <- c(100, 250, 500, 1000, 2500, 5000, 10000)
cum_layers <- calc_zoi_cumulative(x = xx,
                                  radius = rad,
                                  type = "exp_decay",
                                  where = "GRASS",
                                  g_overwrite = TRUE,
                                  verbose = TRUE)
cum_layers <- paste0("hardanger_buff_12km_trails_rast_bin_zoi_cumulative_exp_decay", rad)

xx2 <- sub(pattern = "_bin", replacement = "_1null", xx)
near_layers <- calc_zoi_nearest(x = xx2,
                                radius = rad,
                                type = "exp_decay",
                                where = "GRASS",
                                g_overwrite = TRUE,
                                verbose = TRUE)
near_layers <- paste0("hardanger_buff_12km_trails_rast_1null_zoi_nearest_exp_decay", rad)

rasters_temp <- rgrass::read_RAST(c(cum_layers, near_layers), return_format = "terra")
plot(rasters_temp[[c(2,6,9,13)]])

#--------
# get public roads - skip
# roads_public_all <- sf::st_read(con, DBI::Id(schema = "sam_env", table = "roads_public_renrein_no"))
# # cut
# roads_public <- sf::st_intersection(roads_public_all, reindeer_area) |>
#   dplyr::mutate(value = 1) |>
#   dplyr::select(id, name = gatenavn, publ_priv, traffic_bin, name_area, value)
#
# roads_public
# plot(roads_public[4])
#
# # save
# sf::st_write(roads_public, dsn = "inst/vector/hardanger_roads_public.gpkg", delete_dsn = TRUE)
# # roads_public <- sf::st_read(dsn = "inst/vector/hardanger_roads_public.gpkg")
# # now we document by hand
#
# # try to read
# # cabins count
# (s <- system.file("vector/hardanger_roads_public.gpkg", package = "oneimpact"))
# v <- sf::st_read(s)
# plot(v[4])

#--------
# get private roads - too big file, we skip it
# roads_private_all <- sf::st_read(con, DBI::Id(schema = "sam_env", table = "roads_private_renrein_no"))
# # cut
# roads_private <- sf::st_intersection(roads_private_all, reindeer_area) |>
#   dplyr::mutate(value = 1) |>
#   dplyr::select(id, name = gatenavn, publ_priv, traffic_bin, name_area, value)
#
# roads_private
# plot(roads_private[4])
#
# # save
# sf::st_write(roads_private, dsn = "inst/vector/hardanger_roads_private.gpkg", delete_dsn = TRUE)
# roads_private <- sf::st_write(dsn = "inst/vector/hardanger_roads_private.gpkg")
# # now we document by hand
#
# # try to read
# # cabins count
# (s <- system.file("vector/hardanger_roads_private.gpkg", package = "oneimpact"))
# v <- sf::st_read(s)
# plot(v[4])

#----------------------------
# re-export rasters for cabins private, cabins public, and trails

# region - only hardanger now
rgrass::execGRASS("r.mask", flags = "r")

rgrass::execGRASS("g.region", vector = "hardanger_temp_vector", flags = c("a", "p"))
rgrass::execGRASS("r.mask", vector = "hardanger_temp_vector")

# variables names in GRASS GIS
layers <- rgrass::execGRASS("g.list", type = "raster", pattern = "hardanger_buff_12km*", mapset = ms) |>
  attr("resOut")

# find the layers of interest
priv_cab_g <- layers |>
  grep(pattern = "cabins_private_1null_zoi_nearest|cabins_private_count_zoi_cumulative", value = TRUE) |>
  grep(pattern = "_bin$|count$|1null$", value = TRUE, invert = TRUE)
priv_cab_g <- priv_cab_g[c(1,4,6,2,5,7,3, 7+c(1,4,6,2,5,7,3))]
pub_cab_high_g <- layers |>
  grep(pattern = "cabins_public_hotels_1null_zoi_nearest|cabins_public_hotels_count_zoi_cumulative", value = TRUE) |>
  grep(pattern = "_bin$|count$|1null$", value = TRUE, invert = TRUE)
pub_cab_high_g <- pub_cab_high_g[c(1,4,6,2,5,7,3, 7+c(1,4,6,2,5,7,3))]
trails_g <- layers |>
  grep(pattern = "trails_rast_1null_zoi_nearest|trails_rast_bin_zoi_cumulative", value = TRUE) |>
  grep(pattern = "_bin$|count$|1null$", value = TRUE, invert = TRUE)
trails_g <- trails_g[c(1,4,6,2,5,7,3, 7+c(1,4,6,2,5,7,3))]

all_g <- c(priv_cab_g, pub_cab_high_g, trails_g)
rasters <- rgrass::read_RAST(all_g, return_format = "terra")
names(rasters) <- all_g |>
  sub(pattern = "hardanger_buff_12km_", replacement = "") |>
  sub(pattern = "_rast_bin", replacement = "") |>
  sub(pattern = "_zoi", replacement = "") |>
  sub(pattern = "_hotels", replacement = "") |>
  sub(pattern = "_count", replacement = "")

plot(rasters[[c(t(outer(c(0,14,28), c(2,6), "+")))]])

rgrass::execGRASS("r.mask", flags = "r")

# export
terra::writeRaster(rasters, "data-raw/rast_zois_recomputed_hardanger.tif", gdal = c("COMPRESS=DEFLATE"), overwrite = TRUE)

rasters_500 <- terra::aggregate(rasters, fact = 5, fun = "mean")

# read
rasters_500_orig <- terra::rast("inst/raster/hardanger_rast_predictors_500m.tif")
names(rasters_500_orig)
rasters_500_orig_keep <- rasters_500_orig[[c(1:4,37:39)]]
# plot(rasters_500_orig_keep)

# append with the new ZOIs
rasters_500_updated <- c(rasters_500_orig_keep, rasters_500)
names(rasters_500_updated)

names(rasters_500_updated) <- names(rasters_500_updated) |>
  sub(pattern = "1_null_", replacement = "") |>
  sub(pattern = "1null_", replacement = "") |>
  sub(pattern = "rast_", replacement = "")

# re-write
terra::writeRaster(rasters_500_updated, "data-raw/hardanger_rast_predictors_500m.tif", gdal = c("COMPRESS=DEFLATE"),
                   overwrite = TRUE)
file.copy("data-raw/hardanger_rast_predictors_500m.tif", "inst/raster/", overwrite = TRUE)
file.copy("data-raw/hardanger_rast_predictors_500m.tif.aux.xml", "inst/raster/", overwrite = TRUE)

# rast_predictors$NORUTreclass <- as.numeric(rast_predictors$NORUTreclass)
# rast_df <- terra::as.data.frame(rast_predictors, xy = TRUE, cells = TRUE, na.rm = FALSE)
# saveRDS(rast_df, file = "data-raw/rast_predictors_Hardanger.rds")
# arrow::write_parquet(rast_df, "data-raw/rast_predictors_Hardanger.parquet", compression = "zstd")
# arrow::write_feather(rast_df, "data-raw/rast_predictors_Hardanger.arrow", compression = "zstd")

# test
# put these rasters together with those
(f <- system.file("raster/hardanger_rast_predictors_500m.tif", package = "oneimpact"))
r <- terra::rast(f)
names(r)
plot(r[[20]])

# re-annotate the data for the RSF analysis with these layers

# add linear infrastructure to Hardanger RSF example
library(terra)

# prepared data
reindeer_rsf
# saveRDS(reindeer_rsf, "data-raw/reindeer_rsf_old_only_cabins_pca_landuse_hardanger.rda")

# base data
load("/data/P-Prosjekter/41203800_oneimpact/04_tools/support_oneimpact/cuminf_zoi_GPS_dataset_annotated.rda")
dat
str(dat)

# check
nrow(reindeer_rsf)
nrow(dat)

# subset
dat <- dat |>
  dplyr::select(points_id, x33, y33, herd, use)

# annotate
r1 <- terra::rast("inst/raster/hardanger_rast_predictors_500m.tif")
ext1 <- terra::extract(r1, dat[, c("x33", "y33")])

colnames(reindeer_rsf) %in% colnames(ext1)
colnames(ext1) %in% colnames(reindeer_rsf)

# combine
names(ext1)
ext1[,-1] |> names()

reindeer_rsf <- cbind(dat[, 5, drop = FALSE], ext1[,-1])
names(reindeer_rsf)

# export for the package
usethis::use_data(reindeer_rsf, overwrite = TRUE)
