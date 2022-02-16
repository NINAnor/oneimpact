# Running calc_influence_nearest through R
library(mobsim)
library(terra)
library(dplyr)
library(sf)

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
d <- calc_influence_nearest(pts$rast)
plot(d)

# calculate log_dist (the rest is equal)
log_d <- calc_influence_nearest(pts$rast, type = "log", log_base = 10)
plot(log_d)

# calculate sqrt_dist
sqrt_d <- calc_influence_nearest(pts$rast, type = "sqrt")
plot(sqrt_d)

# calculate exponential decay influence using exp_decay_parms parameter
exp_d1 <- calc_influence_nearest(pts$rast, type = "exp_decay", exp_decay_parms = c(1, 0.001))
plot(exp_d1)

# calculate exponential decay influence using zoi and zoi_decay_threshold (default)
zoi2 <- 1000 # zoi = 1000m
zoi_decay_threshold2 <- 0.05 # here zoi is the distance where the function reaches 0.05
exp_d2 <- calc_influence_nearest(pts$rast, type = "exp_decay", zoi = zoi2,
                                 zoi_decay_threshold = zoi_decay_threshold2)
plot(exp_d2)
# buffer
pts_shp <- pts$pts %>%
  sf::st_as_sf(coords = c(1,2))
# zoi = 1000m
pts_shp %>%
  sf::st_buffer(dist = half_life2*zoi_hl_ratio2) %>%
  sf::st_union() %>%
  plot(add = T, border = "black")
legend("bottomright", legend = c("ZoI"), col = c("black"), lwd=1.1)

# calculate exponential decay influence using half life parameter
# if half_life = 250 m and zoi_hl_ratio = 4, zoi is 1000 m
half_life3 <- 250 # intensity gets down to 1/16 = 0.06 for 4*half_life=1000m
zoi_hl_ratio3 <- 4 # default
exp_d4 <- calc_influence_nearest(pts$rast, type = "exp_decay", half_life = half_life3,
                                 zoi_hl_ratio = zoi_hl_ratio3)
plot(exp_d4)
# buffer
pts_shp <- pts$pts %>%
  sf::st_as_sf(coords = c(1,2))
# half_life = 250m
pts_shp %>%
  sf::st_buffer(dist = half_life2) %>%
  sf::st_union() %>%
  plot(add = T, border = "red")
# zoi = 1000m
pts_shp %>%
  sf::st_buffer(dist = half_life2*zoi_hl_ratio2) %>%
  sf::st_union() %>%
  plot(add = T, border = "black")
legend("bottomright", legend = c("half life", "ZoI"), col = c("red", "black"), lwd=1.1)

# calculate exponential decay influence using zoi parameter and zoi_hl_ratio
zoi4 <- 4000 # intensity gets down to 1/16 = 0.06 for zoi = 4000m, half_life = 1000m
zoi_hl_ratio4 <- 6 # default
exp_d4 <- calc_influence_nearest(pts$rast, type = "exp_decay", zoi = zoi3,
                                 zoi_hl_ratio = zoi_hl_ratio4)
plot(exp_d4)
# buffer
pts_shp <- pts$pts %>%
  sf::st_as_sf(coords = c(1,2))
# half_life = 1000m
pts_shp %>%
  sf::st_buffer(dist = zoi4/zoi_hl_ratio4) %>%
  sf::st_union() %>%
  plot(add = T, border = "red")
# zoi = 4000m
pts_shp %>%
  sf::st_buffer(dist = zoi4) %>%
  sf::st_union() %>%
  plot(add = T, border = "black")
legend("bottomright", legend = c("half life", "ZoI"), col = c("red", "black"), lwd=1.1)

# bartlett influence, ZOI = 2000m
bart_d <- calc_influence_nearest(pts$rast, type = "bartlett", zoi = 2000)
plot(bart_d)

# buffer 2000m
pts_shp %>%
  sf::st_buffer(dist = 2000) %>%
  sf::st_union() %>%
  plot(add = T, border = "black")
legend("bottomright", legend = c("ZoI"), col = c("black"), lwd=1.1)

# calculate threshold influence
d <- calc_influence_nearest(pts$rast, type = "threshold", zoi = 2000)
plot(d)

# Gaussian decay influence
g_d <- calc_influence_nearest(pts$rast, type = "Gauss", zoi = 2000)
plot(g_d)

# buffer 2000m
pts_shp %>%
  sf::st_buffer(dist = 2000) %>%
  sf::st_union() %>%
  plot(add = T, border = "black")
legend("bottomright", legend = c("ZoI"), col = c("black"), lwd=1.1)

#--------------------
