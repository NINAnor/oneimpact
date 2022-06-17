f <- system.file("raster/cabins.tif", package = "oneimpact")
cabins <- terra::rast(f)

# calculate influence
ni <- calc_influence_nearest(cabins, zoi = 1000, type = "exp_decay")
ci <- calc_influence_cumulative(cabins, zoi = 1000, type = "exp_decay")
plot(c(ni, ci))

# rescale
plot(rescale_vals(ci))
plot(rescale_vals(c(ni, ci)))
plot(rescale_vals(c(ni, ci), to = c(0, 100)))
plot(rescale_vals(c(ni, ci), from = c(0, 50)))
