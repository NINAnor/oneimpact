library(terra)
library(dplyr)

# load example raster in metric system
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f) %>% 
  terra::project("EPSG:32631")
# terra::ext(r)[1:2] %>% diff

# create exponential filter
filt_exp1000 <- create_filter(r, zoi = 1000, method = "exp_decay",
                              max_dist = 5000,
                              normalize = T)
filt_exp3000 <- create_filter(r, zoi = 3000, method = "exp_decay", 
                              max_dist = 5000,
                              normalize = T)
# use exponential filter
neigh_r_exp1000 <- terra::focal(r, filt_exp1000, fun = "sum", 
                                na.policy="omit", na.rm=TRUE)
neigh_r_exp3000 <- terra::focal(r, filt_exp3000, fun = "sum", 
                                na.policy="omit", na.rm=TRUE)

# plot
plot(c(r, neigh_r_exp1000, neigh_r_exp3000),
     main = c("original", "exp filter 1000m", "exp filter 3000m"))

# create step filter
filt_step3000 <- create_filter(r, zoi = 3000, method = "step",
                               normalize = T)
# use step filter
neigh_r_step3000 <- terra::focal(r, filt_step3000, fun = "sum", 
                                 na.policy="omit", na.rm=TRUE)

# plot
plot(c(neigh_r_exp3000, neigh_r_step3000), 
     main = c("exp filter 3000m", "step filter 3000m"))
# plot(app(c(neigh_r_exp3000, neigh_r_step3000), "diff"))

# create bartlett filter
filt_bart3000 <- create_filter(r, zoi = 3000, method = "bartlett",
                               normalize = T)
# use bartlett filter
neigh_r_bart3000 <- terra::focal(r, filt_bart3000, fun = "sum",
                                 na.policy="omit", na.rm=TRUE)

# plot
plot(c(neigh_r_exp3000, neigh_r_step3000, neigh_r_bart3000), 
     main = c("exp filter 3000m", "step filter 3000m", "Bartlett filter 3000m"))
# plot(app(c(neigh_r_exp3000, neigh_r_bart3000), "diff"))
# plot(app(c(neigh_r_step3000, neigh_r_bart3000), "diff"))

# save outside R for use in GRASS GIS
# create_filter(r, zoi = 1000, method = "bartlett",
#               max_dist = 5000,
#               normalize = T, save_txt = TRUE)
