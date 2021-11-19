library(mobsim)
library(raster)

set.seed(1234)

# set points
ext <- 30000
wd <- ext/20
pts <- set_points(n_features = 1000, centers = 1,
                  width = wd, res = 100,
                  extent_x = c(0, ext), extent_y = c(0, ext))
plot(pts$pts)
plot(pts$rast)

# calculate distance to the nearest feature
d <- calc_dist(pts$rast)
plot(d)

# calculate log_dist (the rest is equal)
log_d <- calc_dist(pts$rast,
                   transform_dist = "log", log_base = 10)
plot(log_d)

# calculate sqrt_dist (the rest is equal)
sqrt_d <- calc_dist(pts$rast, transform_dist = "sqrt")
plot(sqrt_d)

