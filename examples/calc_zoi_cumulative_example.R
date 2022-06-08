# Running calc_zoi_cumulative through R
library(mobsim)
library(dplyr)
library(terra)

set.seed(1234)

# set points
ext <- 30000
wd <- ext/20
pts <- set_points(n_features = 1000, centers = 1,
                  width = wd, res = 100,
                  extent_x = c(0, ext), extent_y = c(0, ext),
                  buffer_around = 10000)#, use_terra = F)
plot(pts$pts)
plot(pts$rast)

# calculate cumulative zone of influence for multiple influence radii,
# using a Gaussian filter
zoi_values <- c(250, 500, 1000, 2500, 5000)
cumzoi_gauss <- calc_zoi_cumulative(pts$rast, type = "Gauss", zoi_radius = zoi_values,
                                    extent_x_cut = c(0, ext), extent_y_cut = c(0, ext))
plot(cumzoi_gauss)

# calculate cumulative zone of influence for multiple influence radii,
# using a circle neighborhood
cumzoi_circle <- calc_zoi_cumulative(pts$rast, type = "circle", zoi_radius = zoi_values,
                                     extent_x_cut = c(0, ext), extent_y_cut = c(0, ext))
plot(cumzoi_circle)

# calculate cumulative zone of influence for multiple influence radii,
# using an exponential decay neighborhood
cumzoi_exp <- calc_zoi_cumulative(pts$rast, type = "exp_decay", zoi_radius = zoi_values,
                                  extent_x_cut = c(0, ext), extent_y_cut = c(0, ext))
plot(cumzoi_exp)

# comparing
plot(c(cumzoi_gauss[[3]], cumzoi_circle[[3]], cumzoi_exp[[3]]),
     main = c("Gaussian 1000m",
              "Circle 1000m",
              "Exponential decay 1000m"))

# calculate cumulative influence for a single zone of influence
# using a user-defined filter
my_filter <- create_filter(pts$rast, zoi_radius = 1000, type = "rectangle")
cumzoi_user <- calc_zoi_cumulative(pts$rast, type = "mfilter", zoi_radius = my_filter,
                                   extent_x_cut = c(0, ext), extent_y_cut = c(0, ext))
plot(cumzoi_user,
     main = "User-defined rectangular filter")

# calculate density with 1000m radius using an exp_decay neighborhood
density_exp <- calc_zoi_cumulative(pts$rast, type = "exp_decay", zoi_radius = 1000,
                                   output_type = "density",
                                   extent_x_cut = c(0, ext), extent_y_cut = c(0, ext))
# compare
# note the difference in the color scales
plot(c(cumzoi_exp[[3]], density_exp),
     main = c("Cumulative ZoI 1000m", "Density 1000m"))

#--------------------
