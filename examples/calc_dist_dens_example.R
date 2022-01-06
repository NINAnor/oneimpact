library(mobsim)
library(raster)

set.seed(1234)

# set points
ext <- 30000
wd <- ext/20
pts <- set_points(n_features = 1000, centers = 1,
                  width = wd, res = 100,
                  extent_x = c(0, ext), extent_y = c(0, ext),
                  buffer_around = 10000)
plot(pts$pts)
plot(pts$rast)

# calculate distance and densities considering only the initial extent
scales <- c(250, 500, 1000, 2500, 5000)/2
dist_dens <- calc_dist_dens(pts$rast, type_density = "Gauss", scale = scales,
                            extent_x_cut = c(0, ext), extent_y_cut = c(0, ext))
plot(dist_dens)

# calculate exp_decay (the rest is equal)
log_dist_dens <- calc_dist_dens(pts$rast, type_density = "Gauss", scale = scales,
                                transform_dist = "exp_decay", half_life = 200,
                                extent_x_cut = c(0, ext), extent_y_cut = c(0, ext))
plot(log_dist_dens)


# set points
ext <- 30000
wd <- ext/20
pts <- set_points(n_features = 1000, centers = 1,
                  width = wd, res = 100,
                  extent_x = c(0, ext), extent_y = c(0, ext),
                  buffer_around = 0)
plot(pts$pts)
plot(pts$rast)

# calculate distance and densities considering only the initial extent
scales <- c(250, 500, 1000, 2500, 5000)/2
dist_dens <- calc_dist_dens(pts$rast, type_density = "Gauss", scale = scales, pad = T)
plot(dist_dens)
