library(rgrass)
library(terra)
library(arrow)
library(oneimpact)

# connect to GRASS
ms <- "u_bb_cuminf"
NinaR::grassConnect(mapset = ms)

source("/data/P-Prosjekter/41203800_oneimpact/05_papers/02_cumulative_zoi_paper/notebooks/91_find_layers_GRASS.R")

# region
rgrass::execGRASS("g.region", vector = "study_area", flags = c("a", "p"))
# mask
rgrass::execGRASS("r.mask", vector = "study_area", flags = "overwrite")

# rasters - additional ones besides the ones in the MEE paper

# variables names in GRASS GIS
layers <- rgrass::execGRASS("g.list", type = "raster", pattern = "*inf*", mapset = ms) |>
  attr("resOut")

# find the layers of interest
priv_cab_g <- util_find_layer_GRASS(list("private_cabins", "exp_decay"),
                                    layers_grass = layers) |>
  grep(pattern = "_bin_|30000", value = TRUE, invert = TRUE)
pub_cab_high_g <- util_find_layer_GRASS(list("public_cabins_high", "exp_decay"),
                                        layers_grass = layers) |>
  grep(pattern = "_bin_|30000", value = TRUE, invert = TRUE)
priv_road_g <- util_find_layer_GRASS(list("roads_low", "exp_decay"),
                                    layers_grass = layers) |>
  grep(pattern = "30000", value = TRUE, invert = TRUE)
public_road_g <- util_find_layer_GRASS(list("roads_high", "exp_decay"),
                                     layers_grass = layers) |>
  grep(pattern = "30000", value = TRUE, invert = TRUE)
trails_g <- util_find_layer_GRASS(list("trail", "exp_decay"),
                                       layers_grass = layers) |>
  grep(pattern = "30000", value = TRUE, invert = TRUE)


m_vars_g <- c(priv_cab_g, pub_cab_high_g)
m_vars_g <- c(priv_road_g, public_road_g, trails_g)

# mapsets where variables are located
ms_cuminf <- "u_bb_cuminf"
mapsets <- c(ms_cuminf)

raster_names <- paste(m_vars_g, mapsets, sep = "@")

# retrieve rasters
rasters_cabins <- rgrass::read_RAST(raster_names, return_format = "terra")
names(rasters_cabins) <- m_vars_g

# retrieve rasters - linear infrastructure
rasters_linear <- rgrass::read_RAST(raster_names, return_format = "terra")
names(rasters_linear) <- m_vars_g

# remove mask from GRASS GIS
rgrass::execGRASS("r.mask", flags = "r")

# check
plot(rasters_cabins[[1]])
plot(rasters_linear[[2]])

# stack and save rasters
paste0("'", names(rasters_cabins), "'", collapse = ",")
terra::writeRaster(rasters_cabins, "/data/P-Prosjekter/41203800_oneimpact/05_papers/02_cumulative_zoi_paper/data/analysis_GPS/complementary_rasters_to_predict_rsf_cabins.tif")

# paste0("'", names(rasters_linear), "'", collapse = ",")
terra::writeRaster(rasters_linear, "data-raw/rast_predictors_hardanger_100_linear.tif", gdal = c("COMPRESS=DEFLATE"))

# preparing and saving it
# model
bag_object$formula

# coefficients
m_vars <- all.vars(bag_object$formula)[-1]

# load vectors
path <- "/data/P-Prosjekter/41203800_oneimpact/05_papers/02_cumulative_zoi_paper/data/analysis_GPS/"
study_area_v <- terra::vect(paste0(path, "study_area_v.gpkg"))
private_cabins_v <- terra::vect(paste0(path, "private_cabins_v.gpkg"))
public_cabins_high_v <- terra::vect(paste0(path, "public_cabins_high_v.gpkg"))
use_v <- terra::vect(paste0(path, "use_v.gpkg"))

# load raster
all_rast <- terra::rast("/data/P-Prosjekter/41203800_oneimpact/05_papers/02_cumulative_zoi_paper/data/analysis_GPS/rasters_to_predict_rsf.tif")
names(all_rast) <- c('private_cabins_cumulative_threshold_10000','public_cabins_high_cumulative_exp_decay_20000','norway_pca_klima_axis1','norway_pca_klima_axis2','norway_pca_klima_axis3','norway_pca_klima_axis4','NORUTreclass','private_cabins_inf_nearest_threshold10000','public_cabins_high_inf_nearest_exp_decay20000')
rasters <- all_rast[[3:6]]
lu_rast <- all_rast[[7]]

# cabins
rasters_cabins <- terra::rast("/data/P-Prosjekter/41203800_oneimpact/05_papers/02_cumulative_zoi_paper/data/analysis_GPS/complementary_rasters_to_predict_rsf_cabins.tif")
radii <- c(100, 1000, 10000, 20000, 250, 2500, 500, 5000)
names(rasters_cabins) <- c(paste0("private_cabins_", "cumulative", "_exp_decay_", radii),
                           paste0("private_cabins_", "nearest", "_exp_decay_", radii),
                           paste0("public_cabins_high_", "cumulative", "_exp_decay_", radii),
                           paste0("public_cabins_high_", "nearest", "_exp_decay_", radii))

# prepare rasters for prediction
# scale continuous variables
# rasters_scaled[[1:2]] <- scale(rasters_scaled[[1:2]])
# terra::values(rasters_scaled[[1:2]]) <- scale(terra::values(rasters_scaled[[1:2]]))
# reclassify land cover raster
classes <- data.frame(id = sort(unique(values(lu_rast))),
                      NORUTreclass = as.character(sort(unique(values(lu_rast)))))
classes$NORUTreclass <- ifelse(classes$id < 9, "11forest",
                               ifelse(classes$id < 12, "bog",
                                      ifelse(classes$id < 21, "mountain",
                                             ifelse(classes$id < 22, "glacier",
                                                    ifelse(classes$id < 23, "water",
                                                           ifelse(classes$id < 25,"other", NA))))))
classes$NORUTreclass <- ifelse(classes$NORUTreclass == "mountain", as.character(classes$id), classes$NORUTreclass)
levels(lu_rast) <- as.data.frame(classes)
activeCat(lu_rast)
# check
# sort(unique(values(lu_rast)))
# levels(lu_rast)[[1]]
# sort(unique(dat$NORUTreclass))
plot(lu_rast)
plot(rasters)
plot(rasters_cabins[[1:8]])
plot(rasters_cabins[[9:16]])
plot(rasters_cabins[[17:24]])
plot(rasters_cabins[[25:32]])

# rasters for prediction
rast_predictors <- c(rasters, lu_rast, rasters_cabins[[c(1:8, 17:24)]])

rast_predictors$norway_pca_klima_axis1_sq <- rast_predictors$norway_pca_klima_axis1**2
rast_predictors$norway_pca_klima_axis2_sq <- rast_predictors$norway_pca_klima_axis2**2

# degrade resolution
rast_predictors_500 <- terra::aggregate(rast_predictors[[-5]], fact = 5, fun = "mean")
rast_predictors_500$NORUTreclass <- terra::aggregate(rast_predictors$NORUTreclass, fact = 5, fun = "modal")
plot(rast_predictors_500$NORUTreclass)
rast_predictors_500$NORUTreclass |> levels()

terra::writeRaster(rast_predictors, "data-raw/rast_predictors_hardanger_100.tif", gdal = c("COMPRESS=DEFLATE"))
terra::writeRaster(rast_predictors_500, "inst/raster/rast_predictors_hardanger_500.tif", gdal = c("COMPRESS=DEFLATE"),
                   overwrite = TRUE)
file.copy("data-raw/rast_predictors_hardanger_500.tif", "inst/raster/", overwrite = TRUE)

rast_predictors$NORUTreclass <- as.numeric(rast_predictors$NORUTreclass)
rast_df <- terra::as.data.frame(rast_predictors, xy = TRUE, cells = TRUE, na.rm = FALSE)
saveRDS(rast_df, file = "data-raw/rast_predictors_Hardanger.rds")
arrow::write_parquet(rast_df, "data-raw/rast_predictors_Hardanger.parquet", compression = "zstd")
arrow::write_feather(rast_df, "data-raw/rast_predictors_Hardanger.arrow", compression = "zstd")

# test
rast_df <- readRDS("data-raw/rast_predictors_Hardanger.rds")
str(rast_df)
rast <- terra::rast(rast_df[,-1], type = "xyz", crs = "epsg:25833")
rast

plot(rast)

rast_predictors <- terra::rast("data-raw/rast_predictors_hardanger_100.tif")

