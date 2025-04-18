#-----
# using mobsim
library(terra)
library(mobsim)

set.seed(1234)

# gradient distribution
ext <- 30000
wd <- ext/5
pts <- set_points(n_features = 1000, centers = 1,
                  width = wd, res = 100,
                  extent_x = c(0, ext), extent_y = c(0, ext))
plot(pts$pts)
plot(pts$rast, col = "black")

# one focus of features, with buffer around
wd <- ext/20
pts <- set_points(n_features = 1000, centers = 1,
                  width = wd, res = 100,
                  extent_x = c(0, ext), extent_y = c(0, ext),
                  buffer_around = 10000)
plot(pts$pts)
plot(pts$rast, col = "black")

#-----
# using base raster

# raster
set.seed(12)
r <- raster::raster(matrix(runif(12),3,4)) |>
  raster::disaggregate(fact = 10)

# points from raster
pts <- set_points(n_features = 100, method = "raster",
                  base_raster = r)
plot(pts$base_rast)
plot(pts$pts)
plot(pts$rast, col = "black")

#-----
# using random or regular

set.seed(123)
ext <- 30000
pts <- set_points(n_features = 1000, method = "random",
                  res = 100,
                  extent_x = c(0, ext), extent_y = c(0, ext))
plot(pts$pts)
plot(pts$rast, col = "black")

pts <- set_points(n_features = 1000, method = "regular",
                  res = 100,
                  extent_x = c(0, ext), extent_y = c(0, ext))
plot(pts$pts)
plot(pts$rast, col = "black")

#-----
# using point coordinates as input
pt_input <- data.frame(x = c(0.5, 0.7), y = c(0.5, 0.3))
pts <- set_points(point_coordinates = pt_input)
plot(pts$pts)
plot(pts$rast, col = "black")
