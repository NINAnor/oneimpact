#-----
# minimal example

# example based on
# https://gis.stackexchange.com/questions/224321/randomly-generate-points-using-weights-from-raster
library(raster)

# raster
set.seed(12)
r <- raster::raster(matrix(runif(12),3,4))

# points
pts <- set_points_from_raster(r, n_features = 300)

# plot
raster::plot(r)
points(pts)

# with terra
r <- terra::rast(r)
# points
pts <- set_points_from_raster(r, n_features = 300)
