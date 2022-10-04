# Running calc_zoi_cumulative through R
library(terra)

# Load raster data
f <- system.file("raster/sample_area_cabins_count.tif", package = "oneimpact")
cabins <- terra::rast(f)

# calculate cumulative zone of influence for multiple influence radii,
# using a Gaussian filter
zoi_values <- c(250, 500, 1000, 2500, 5000)
cumzoi_gauss <- calc_zoi_cumulative(cabins, type = "Gauss", radius = zoi_values)
plot(cumzoi_gauss)

# calculate cumulative zone of influence for multiple influence radii,
# using a circle neighborhood
cumzoi_circle <- calc_zoi_cumulative(cabins, type = "circle", radius = zoi_values)
plot(cumzoi_circle)

# calculate cumulative zone of influence for multiple influence radii,
# using an exponential decay neighborhood
cumzoi_exp <- calc_zoi_cumulative(cabins, type = "exp_decay", radius = zoi_values)
plot(cumzoi_exp)

# comparing
plot(c(cumzoi_gauss[[3]], cumzoi_circle[[3]], cumzoi_exp[[3]]),
     main = c("Gaussian 1000m",
              "Circle 1000m",
              "Exponential decay 1000m"))

# calculate cumulative influence for a single zone of influence
# using a user-defined filter
my_filter <- filter_create(cabins, radius = 1000, type = "rectangle")
cumzoi_user <- calc_zoi_cumulative(cabins, type = "mfilter", radius = my_filter)
plot(cumzoi_user, main = "User-defined rectangular filter")

# calculate density with 1000m radius using an exp_decay neighborhood
density_exp <- calc_zoi_cumulative(cabins, type = "exp_decay", radius = 1000,
                                   output_type = "density")
# compare
# note the difference in the color scales
plot(c(cumzoi_exp[[3]], density_exp),
     main = c("Cumulative ZoI 1000m", "Density 1000m"))

#--------------------
