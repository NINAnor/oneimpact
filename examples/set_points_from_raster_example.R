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

raster::plot(r)
points(pts)

#-----
# using NLMR
library(NLMR)

# example NLM
set.seed(123)
nlm1 <- NLMR::nlm_mpd(100, 100, 100, roughness = .5)

# points
pts <- set_points_from_raster(nlm1, n_features = 1000)

# plot
raster::plot(nlm1)
points(pts)

