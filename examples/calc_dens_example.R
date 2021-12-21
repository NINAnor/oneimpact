library(mobsim)
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

# calculate densities for multiple scales, 
# considering only the initial extent, for a Gaussian filter
scales <- c(250, 500, 1000, 2500, 5000)/2
dens <- calc_dens(pts$rast, type = "Gauss", scale = scales,
                  extent_x_cut = c(0, ext), extent_y_cut = c(0, ext))
plot(dens)

# calculate densities for multiple scales, 
# using a circle neighborhood
dens_circle <- calc_dens(pts$rast, type = "circle", scale = scales,
                         extent_x_cut = c(0, ext), extent_y_cut = c(0, ext))
plot(dens_circle)

# calculate densities for multiple scales, 
# using a user-defined filter (exp_decay)
my_filter <- create_filter(zoi = 1000, res = 100, method = "exp_decay")
dens_exp <- calc_dens(pts$rast, type = "mfilter", scale = my_filter,
                         extent_x_cut = c(0, ext), extent_y_cut = c(0, ext))
plot(c(dens[[3]], dens_circle[[3]], dens_exp),
     main = c("Gaussian filter", "Circle neighborhood", "Exponential decay filter"))

     