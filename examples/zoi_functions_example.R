# generic dist_decay function
oneimpact::dist_decay(500, radius = 1000, type = "exp_decay")
oneimpact::dist_decay(500, radius = 1000, type = "gaussian_decay")
oneimpact::dist_decay(500, radius = 1000, type = "linear_decay")
oneimpact::dist_decay(500, radius = 1000, type = "step_decay")

# test the zone of influence functions
# here we use ggplot() to illustrate the functions, to make the figures more
# widely reproducible
# to ease the plots, use the function plot_zoi1d()
library(ggplot2)

# exponential decay
exp_decay(10, radius = 30)

f1 <- ggplot(data.frame(x = c(0, 30)), aes(x = x)) +
  stat_function(fun = exp_decay, args = list(radius = 20)) +
  labs(x = "Distance", y = "Zone of Influence") +
  theme_bw()
f1

# exponential decay - two sides
f1_2 <- ggplot(data.frame(x = c(-30, 30)), aes(x = x)) +
  stat_function(fun = exp_decay,
                args = list(radius = 20, oneside = FALSE)) +
  labs(x = "Distance", y = "Zone of Influence") +
  theme_bw()
f1_2

# threshold
threshold_decay(5, radius = 10)
threshold_decay(10, radius = 10)

f2 <- ggplot(data.frame(x = c(0, 30)), aes(x = x)) +
  stat_function(fun = threshold_decay,
                args = list(radius = 20), linetype = 2) +
  labs(x = "Distance", y = "Zone of Influence") +
  theme_bw()
f2

# threshold - two sides
f2_2 <- ggplot(data.frame(x = c(-30, 50)), aes(x = x)) +
  stat_function(fun = threshold_decay,
                args = list(radius = 20, oneside = FALSE), linetype = 2) +
  labs(x = "Distance", y = "Zone of Influence") +
  theme_bw()
f2_2

# linear, tent, or bartlett decay
bartlett_decay(5, radius = 10)
bartlett_decay(8, radius = 10)

f3 <- ggplot(data.frame(x = c(0, 30)), aes(x = x)) +
  stat_function(fun = bartlett_decay, args = list(radius = 20), linetype = 3) +
  labs(x = "Distance", y = "Zone of Influence") +
  theme_bw()
f3

# linear, two sides
f3_3 <- ggplot(data.frame(x = c(-30, 40)), aes(x = x)) +
  stat_function(fun = bartlett_decay,
                args = list(radius = 20, origin = 10, oneside = FALSE), linetype = 3) +
  labs(x = "Distance", y = "Zone of Influence") +
  theme_bw()
f3_3

# guassian or half normal
gaussian_decay(5, sigma = 6)

f4 <- ggplot(data.frame(x = c(0, 30)), aes(x = x)) +
  stat_function(fun = gaussian_decay,
                args = list(radius = 20, zoi_limit = 0.05), linetype = 4) +
  labs(x = "Distance", y = "Zone of Influence") +
  geom_vline(xintercept = 20, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0.05, linetype = 2, color = "darkgrey") +
  theme_bw()
f4

# half normal - two sides
gaussian_decay(5, sigma = 6)

f4_2 <- ggplot(data.frame(x = c(-30, 30)), aes(x = x)) +
  stat_function(fun = gaussian_decay,
                args = list(radius = 20, zoi_limit = 0.05), linetype = 4) +
  labs(x = "Distance", y = "Zone of Influence") +
  geom_vline(xintercept = c(-20, 20), linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0.05, linetype = 2, color = "darkgrey") +
  theme_bw()
f4_2

# plot several ZoI with the same radius
f1 +
  stat_function(fun = threshold_decay, args = list(radius = 20), linetype = 2) +
  stat_function(fun = bartlett_decay, args = list(radius = 20), linetype = 3) +
  stat_function(fun = gaussian_decay, args = list(radius = 20, zoi_limit = 0.05), linetype = 4) +
  labs(x = "Distance", y = "Zone of Influence") +
  theme_bw()

#---
# applying dist_decay functions for rasters
library(terra)

# calculate Euclidean distance
f <- system.file("raster/sample_area_cabins.tif", package = "oneimpact")
cabins <- terra::rast(f)
cabins_dist <- calc_zoi_nearest(cabins, type = "euclidean")

# transform Euclidean in distance decay
# exponential decay
plot(oneimpact::dist_decay(cabins_dist, radius = 1000, type = "exp_decay"))
# linear decay
plot(oneimpact::dist_decay(cabins_dist, radius = 1000, type = "tent_decay"))
