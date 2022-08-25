library(terra)

# load example - raster of tourist private cabins
f <- system.file("raster/cabins.tif", package="oneimpact")
r <- rast(f)
# terra::ext(r)[1:2] %>% diff

# set value zero where there are no cabins
r[is.na(r)] <- 0

# create exponential filter
filt_exp1000 <- filter_create(r, radius = 1000,
                              zoi_limit = 0.01,
                              type = "exp_decay",
                              max_dist = 5000,
                              normalize = T)
filt_exp3000 <- filter_create(r, radius = 3000,
                              zoi_limit = 0.01,
                              type = "exp_decay",
                              max_dist = 5000,
                              normalize = T)
# use exponential filter
neigh_r_exp1000 <- terra::focal(r, filt_exp1000, fun = "sum",
                                na.policy = "omit", na.rm = TRUE)
neigh_r_exp3000 <- terra::focal(r, filt_exp3000, fun = "sum",
                                na.policy = "omit", na.rm = TRUE)

# plot
plot(c(r, neigh_r_exp1000, neigh_r_exp3000),
     main = c("original", "exp filter 1000m", "exp filter 3000m"))

# create step filter
filt_step3000 <- filter_create(r, radius = 3000, type = "step",
                               normalize = T)
# use step filter
neigh_r_step3000 <- terra::focal(r, filt_step3000, fun = "sum",
                                 na.policy = "omit", na.rm = TRUE)

# plot
plot(c(neigh_r_exp3000, neigh_r_step3000),
     main = c("exp filter 3000m", "step filter 3000m"))
# plot(app(c(neigh_r_exp3000, neigh_r_step3000), "diff"))

# create bartlett (linear/tent decay) filter
filt_bart3000 <- filter_create(r, radius = 3000, type = "bartlett",
                               normalize = T)
# use bartlett filter
neigh_r_bart3000 <- terra::focal(r, filt_bart3000, fun = "sum",
                                 na.policy = "omit", na.rm = TRUE)

# create Gaussian filter - parameterized with zoi
filt_gauss3000 <- filter_create(r, radius = 3000,
                                type = "Gauss",
                                zoi_limit = 0.01,
                                normalize = T)
# use Gaussian filter
neigh_r_gauss3000 <- terra::focal(r, filt_gauss3000, fun = "sum",
                                  na.policy = "omit", na.rm = TRUE)

# plot
plot(c(neigh_r_exp3000, neigh_r_step3000, neigh_r_bart3000, neigh_r_gauss3000),
     main = c("exp filter 3000m", "step filter 3000m",
              "Bartlett filter 3000m", "Gaussian filter 3000m"))
# plot(app(c(neigh_r_exp3000, neigh_r_bart3000), "diff"))
# plot(app(c(neigh_r_step3000, neigh_r_bart3000), "diff"))

# Not run
# save outside R for use in GRASS GIS
if(FALSE) {
  filter_create(r, radius = 1000,
                type = "bartlett",
                max_dist = 5000,
                normalize = T, save_txt = TRUE)
}
