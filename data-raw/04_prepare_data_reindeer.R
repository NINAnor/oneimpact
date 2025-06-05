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
