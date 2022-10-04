library(terra)

f <- system.file("raster/sample_area_cabins.tif", package = "oneimpact")
cabins <- terra::rast(f)

# calculate influence
ni <- calc_zoi_nearest(cabins, radius = 1000, type = "exp_decay")
ci <- calc_zoi_cumulative(cabins, radius = 1000, type = "exp_decay",
                          zeroAsNA = TRUE)
plot(c(ni, ci))

# rescale
plot(raster_rescale(ci)) # rescale to [0,1]
plot(raster_rescale(c(ni, ci))) # rescale both to [0,1]
plot(raster_rescale(c(ni, ci), to = c(0, 100))) # rescale to [0,100]
plot(raster_rescale(c(ni, ci), from = c(0, 50))) # rescale to [0,1] from [0,50]
