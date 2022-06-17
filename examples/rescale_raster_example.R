f <- system.file("raster/cabins.tif", package = "oneimpact")
cabins <- terra::rast(f)

# calculate influence
ni <- calc_zoi_nearest(cabins, zoi_radius = 1000, type = "exp_decay")
ci <- calc_zoi_cumulative(cabins, zoi_radius = 1000, type = "exp_decay",
                          zeroAsNA = TRUE)
plot(c(ni, ci))

# rescale
plot(rescale_raster(ci)) # rescale to [0,1]
plot(rescale_raster(c(ni, ci))) # rescale both to [0,1]
plot(rescale_raster(c(ni, ci), to = c(0, 100))) # rescale to [0,100]
plot(rescale_raster(c(ni, ci), from = c(0, 50))) # rescale to [0,1] from [0,50]
