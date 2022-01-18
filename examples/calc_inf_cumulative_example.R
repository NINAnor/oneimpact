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

# calculate cumulative influence for multiple zones of influence, 
# considering only the initial extent, for a Gaussian filter
zoi_values <- c(250, 500, 1000, 2500, 5000)/2
cuminf <- calc_influence_cumulative(pts$rast, type = "Gauss", zoi = zoi_values,
                                    extent_x_cut = c(0, ext), extent_y_cut = c(0, ext))
plot(cuminf)

# calculate cumulative influence for multiple zones of influence, 
# using a circle neighborhood
zoi_values <- c(250, 500, 1000, 2500, 5000)
cuminf_circle <- calc_influence_cumulative(pts$rast, type = "circle", zoi = zoi_values,
                         extent_x_cut = c(0, ext), extent_y_cut = c(0, ext))
plot(cuminf_circle)

# calculate cumulative influence for multiple zones of influence, 
# using a user-defined filter (exp_decay)
my_filter <- create_filter(zoi = 1000, res = 100, method = "exp_decay")
cuminf_exp <- calc_influence_cumulative(pts$rast, type = "mfilter", zoi = my_filter,
                                        extent_x_cut = c(0, ext), extent_y_cut = c(0, ext))
plot(c(cuminf[[3]], cuminf_circle[[3]], cuminf_exp),
     main = c("Gaussian filter", "Circle neighborhood", "Exponential decay filter"))

     