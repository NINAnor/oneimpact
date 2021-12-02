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

# calculate sqrt_dist
sqrt_d <- calc_dist(pts$rast, transform_dist = "sqrt")
plot(sqrt_d)

# calculate exponential decay distance
exp_d1 <- calc_dist(pts$rast, transform_dist = "exp_decay", exp_decay_parms = c(1, 0.001))
plot(exp_d1)
# using half life
half_life <- 1000 # intensity gets down to 1/16 = 0.06 for 4*half_life
exp_d2 <- calc_dist(pts$rast, transform_dist = "exp_decay", exp_hl = 1000)
plot(exp_d2)
# buffer
pts_shp <- pts$pts %>% 
  sf::st_as_sf(coords = c(1,2))
# 1000m
pts_shp %>% 
  sf::st_buffer(dist = 1000) %>% 
  sf::st_union() %>% 
  plot(add = T)
# 4000m
pts_shp %>% 
  sf::st_buffer(dist = 2000) %>% 
  sf::st_union() %>% 
  plot(add = T)

