library(terra)

# Load raster data
f <- system.file("raster/sample_area_cabins.tif", package = "oneimpact")
cabins <- terra::rast(f)

#---
# check background values
terra::freq(cabins) ## No zeros, background is NA

# compute both Zoi metrics with Gaussian decay, radius = 1000 m
# since the background is NA, we use zeroAsNA = FALSE
zoi_metrics <- calc_zoi(cabins,
                        radius = 1000,
                        type = "Gauss",
                        zeroAsNA = FALSE)
# check
zoi_metrics
# plot
plot(zoi_metrics)

#-------

# Load raster data
f <- system.file("raster/sample_area_cabins_count.tif", package = "oneimpact")
cabins_count <- terra::rast(f)

# check background values
terra::freq(cabins_count) ## Places with no infrastructure have value zero

# compute both Zoi metrics with linear decay, varying radius from 1000 m to 3000 m
# since the background is zero, we use zeroAsNA = TRUE
zoi_metrics2 <- calc_zoi(cabins_count,
                         radius = c(1000, 2000, 3000),
                         type = "bartlett",
                         zeroAsNA = TRUE,
                         output_type = "density")

# check
zoi_metrics2
# plot
plot(zoi_metrics2)
