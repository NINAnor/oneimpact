library(mobsim)
library(raster)

set.seed(1234)

# gradient distribution
ext <- 30000
wd <- ext/5
pts <- set_points(n_features = 1000, centers = 1,
                  width = wd, res = 100,
                  extent_x = c(0, ext), extent_y = c(0, ext))
plot(pts$pts)
plot(pts$rast)

# one focus of features, with buffer around
wd <- ext/20
pts <- set_points(n_features = 1000, centers = 1,
                  width = wd, res = 100,
                  extent_x = c(0, ext), extent_y = c(0, ext),
                  buffer_around = 10000)
plot(pts$pts)
plot(pts$rast)
