# Running calc_zoi_nearest through R
library(terra)
library(sf)

# Load raster data
f <- system.file("raster/sample_area_cabins.tif", package = "oneimpact")
cabins <- terra::rast(f)

# Load vector data
f2 <- system.file("vector/sample_area_cabins.gpkg", package = "oneimpact")
cabins_vect <- sf::st_read(f2)

# calculate distance to the nearest feature
d <- calc_zoi_nearest(cabins, type = "euclidean")
plot(d)

# calculate log_dist (the rest is equal)
log_d <- calc_zoi_nearest(cabins, type = "log", log_base = 10)
plot(log_d)

# calculate sqrt_dist
sqrt_d <- calc_zoi_nearest(cabins, type = "sqrt")
plot(sqrt_d)

# calculate exponential decay zone of influence
# using exp_decay_parms parameter
exp_d1 <- calc_zoi_nearest(cabins, type = "exp_decay",
                           intercept = 1, lambda = 0.001)
plot(exp_d1)

# calculate exponential decay zone of influence using
# radius and zoi_limit (default)
radius2 <- 1000 # zoi = 1000m
zoi_limit2 <- 0.05 # here zoi is the distance where the function reaches 0.05
exp_d2 <- calc_zoi_nearest(cabins, type = "exp_decay", radius = radius2,
                           zoi_limit = zoi_limit2)
plot(exp_d2)
# buffer
# zoi = 1000m
cabins_vect |>
  sf::st_buffer(dist = radius2) |>
  sf::st_union() |>
  plot(add = T, border = "black")
legend("bottomright", legend = c("ZoI radius"), col = c("black"), lwd = 1.1)

# calculate exponential decay zone of influence using half life parameter
# if half_life = 250 m and zoi_hl_ratio = 4, zoi is 1000 m
half_life3 <- 250 # intensity gets down to 1/16 = 0.06 for 4*half_life=1000m
zoi_hl_ratio3 <- 4 # default
exp_d4 <- calc_zoi_nearest(cabins, type = "exp_decay", half_life = half_life3,
                           zoi_hl_ratio = zoi_hl_ratio3)
plot(exp_d4)
# buffer
cabins_vect |>
  sf::st_buffer(dist = half_life3) |>
  sf::st_union() |>
  plot(add = T, border = "red")
# zoi = 1000m
cabins_vect |>
  sf::st_buffer(dist = half_life3*zoi_hl_ratio3) |>
  sf::st_union() |>
  plot(add = T, border = "black")
legend("bottomright", legend = c("Exponential half-life", "ZoI radius"),
       col = c("red", "black"), lwd = 1.1)

# calculate exponential decay zone of influence using
# radius parameter and zoi_hl_ratio
radius4 <- 4000 # intensity gets down to 1/16 = 0.06 for zoi = 4000m, half_life = 1000m
zoi_hl_ratio4 <- 6 # default
exp_d4 <- calc_zoi_nearest(cabins, type = "exp_decay", radius = radius4,
                           zoi_hl_ratio = zoi_hl_ratio4)
plot(exp_d4)
# buffer
# half_life = 1000m
cabins_vect |>
  sf::st_buffer(dist = radius4/zoi_hl_ratio4) |>
  sf::st_union() |>
  plot(add = T, border = "red")
# zoi = 4000m
cabins_vect |>
  sf::st_buffer(dist = radius4) |>
  sf::st_union() |>
  plot(add = T, border = "black", )
legend("bottomright", legend = c("Exponential half-life", "ZoI radius"),
       col = c("red", "black"), lwd = 1.1)

#---
# bartlett influence, ZOI = 1000m
bart_d <- calc_zoi_nearest(cabins, type = "bartlett", radius = 1000)
plot(bart_d)

# buffer 1000m
cabins_vect |>
  sf::st_buffer(dist = 1000) |>
  sf::st_union() |>
  plot(add = T, border = "black")
legend("bottomright", legend = c("Bartlett ZoI 1000m"),
       col = c("black"), lwd = 1.1)

# calculate threshold influence, ZoI = 1000m
d <- calc_zoi_nearest(cabins, type = "threshold", radius = 1000)
plot(d)

# Gaussian decay influence, ZoI = 1000m
g_d <- calc_zoi_nearest(cabins, type = "Gauss", radius = 1000)
plot(g_d)

# buffer 1000m
cabins_vect |>
  sf::st_buffer(dist = 1000) |>
  sf::st_union() |>
  plot(add = T, border = "black")
legend("bottomright", legend = c("Gaussian ZoI 1000m"),
       col = c("black"), lwd = 1.1)

#--------------------
