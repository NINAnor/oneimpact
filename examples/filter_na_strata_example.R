data(reindeer)
# use-availability setup for ssf
r <- reindeer |>
  dplyr::rename(x_ = x, y_ = y, t_ = t) |>
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
  terra::buffer(max(r$sl_)*1.1)
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
plot(cabins_zoi)
# roads
reindeer_roads_v <- terra::vect(system.file("vector/reindeer_roads_public.gpkg", package = "oneimpact"))
reindeer_roads <- terra::rasterize(reindeer_roads_v, rr, field = "value")
roads_zoi <- calc_zoi_cumulative(reindeer_roads, radius = c(seq(1000, 3000, 1000)),
                                 type = "Gauss", zeroAsNA = TRUE)
names(roads_zoi) <- paste0("roads", seq(1000, 3000, 1000))
plot(roads_zoi)

# extract
st <- terra::extract(c(cabins_zoi, roads_zoi), terra::vect(r, geom = c("x1_", "y1_")))
names(st)[-1] <- paste0("startpt_", names(st)[-1])
end <- terra::extract(c(cabins_zoi, roads_zoi), terra::vect(r, geom = c("x2_", "y2_")))
names(end)[-1] <- paste0("endpt_", names(end)[-1])

data_annotated <- dplyr::bind_cols(r, st[-1], end[-1])
data_annotated[grep("cabins", names(data_annotated))[1]:ncol(data_annotated)] <-
  apply(data_annotated[grep("cabins", names(data_annotated))[1]:ncol(data_annotated)], 2, scale)
data_annotated$case_ <- ifelse(data_annotated$case_, 1, 0)

f <- case_ ~ strata(step_id) + sl_*startpt_roadsXXX + sl_*startpt_cabinsXXX +
  endpt_roadsXXX + endpt_cabinsXXX
f <- add_zoi_formula(f, zoi_radius = c(1000, 2000, 3000), pattern = "XXX")
filter_na_strata(f, data_annotated)
# ?amt::remove_incomplete_strata()

data_annotated_v <- terra::vect(data_annotated, geom = c("x2_", "y2_"))
plot(data_annotated_v[data_annotated_v$case_ == F], col = "grey", cex = 0.3)
points(data_annotated_v[data_annotated_v$case_ == T], col = "red", cex = 0.3)
