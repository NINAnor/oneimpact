library(mobsim)
library(raster)

set.seed(1234)

# multiple focii
ext <- 30000
wd <- 0.1*ext
scales <- c(250, 500, 1000, 2500, 5000)/2
rasts <- simulate_dist_dens(n_features = 1000, centers = 5,
                            width = wd, res = 100,
                            extent_x = c(0, ext), extent_y = c(0, ext),
                            buffer_around = 10000,
                            type_density = "Gauss", scale = scales,
                            extent_x_cut = c(0, ext), extent_y_cut = c(0, ext))
plot(rasts)
