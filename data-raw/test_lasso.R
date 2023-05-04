# Libraries ---------------------------
library(terra)
library(oneimpact)

# Load annotated data and prepare for iSSF --------------------

# All pops, steps
filename = c("wrein_ssfdat_sum_20220513.rda")[1]
load(paste0("/data/P-Prosjekter/41203800_oneimpact/04_tools/support_oneimpact/", filename))

str(data)
names(data)

# One pop, pts
load("/data/P-Prosjekter/41203800_oneimpact/04_tools/support_oneimpact/cuminf_zoi_GPS_dataset_annotated.rda")

str(dat)
names(dat)

# Spatial stratification --------------------
# spStrat <- spat_strat(dat[dat$use == 1,], coords = c("x33", "y33"), colH0 = ""),
spStrat <- spat_strat(dat[dat$use == 1,], coords = c("x33", "y33"))

spStrat <- spat_strat(data[data$case == 1,], coords = c("x", "y"), colH0 = "population",
                      colID = "gps_data_animals_id")

spStrat <- spat_strat(data_annotated[data_annotated$case_ == 1,], coords = c("x2_", "y2_"),
                      colH0 = "animal_year_id", colID = "step_id")


# define formula ----------------------------
# f <- case ~ strata(points_id) + step_length + I(step_length^2) +
#   step_length:strtpt_houses_5000 + step_length:strtpt_private_cabins_2500 + step_length:strtpt_pub_cabins_summer_high_2500 + step_length:strtpt_pub_cabins_summer_low_1000 +
#   step_length:strtpt_roads_summer_high_5000 + step_length:strtpt_roads_summer_low_2500 + step_length:strtpt_trails_pseudotui_2500 +  step_length:strtpt_powerlines_1000 + step_length:strtpt_railway_1000 +
#   step_length:strtpt_dem_tpi_500_50m + step_length:strtpt_solar_radiation_10m_july +
#   step_length:strtpt_reservoirs_250 + step_length:strtpt_lakes_250 + step_length:strtpt_spraverageallyrs_500_100m_2 + step_length:I(strtpt_spraverageallyrs_500_100m_2^2) +
#   step_length:strtpt_digestible_biomass_summer+
#   roads_summer_low_100 + roads_summer_high_100 + trails_pseudotui_100 + powerlines_100 + railway_100 +
#   houses_500 + private_cabins_500 + pub_cabins_summer_high_500 + pub_cabins_summer_low_500 +
#   max_slope + reservoirs_250 + lakes_nores_250 + length_water +
#   endpt_houses_5000 + endpt_private_cabins_2500 + endpt_pub_cabins_summer_high_2500 + endpt_pub_cabins_summer_low_1000 +
#   endpt_roads_summer_high_5000 + endpt_roads_summer_low_2500 + endpt_trails_pseudotui_2500 +  endpt_powerlines_1000 + endpt_railway_1000 +
#   endpt_dem_tpi_500_50m + endpt_solar_radiation_10m_july +
#   endpt_reservoirs_250 + endpt_lakes_250 + endpt_spraverageallyrs_500_100m_2 + I(endpt_spraverageallyrs_500_100m_2^2) +
#   endpt_digestible_biomass_summer

f <- case ~ strata(points_id) + step_length + I(step_length^2) +
  step_length:strtpt_houses_5000 + step_length:strtpt_private_cabins_2500 + step_length:strtpt_pub_cabins_summer_high_2500 + step_length:strtpt_pub_cabins_summer_low_1000 +
  step_length:strtpt_roads_summer_high_5000 + step_length:strtpt_roads_summer_low_2500 + step_length:strtpt_trails_pseudotui_2500 +  step_length:strtpt_powerlines_1000 + step_length:strtpt_railway_1000 +
  step_length:strtpt_spraverageallyrs_500_100m_2 + step_length:I(strtpt_spraverageallyrs_500_100m_2^2) +
  roads_summer_low_100 + roads_summer_high_100 + trails_pseudotui_100 + powerlines_100 + railway_100 +
  houses_500 + private_cabins_500 + pub_cabins_summer_high_500 + pub_cabins_summer_low_500 +
  max_slope + reservoirs_250 + length_water +
  endpt_houses_5000 + endpt_private_cabins_2500 + endpt_pub_cabins_summer_high_2500 + endpt_pub_cabins_summer_low_1000 +
  endpt_roads_summer_high_5000 + endpt_roads_summer_low_2500 + endpt_trails_pseudotui_2500 +  endpt_powerlines_1000 + endpt_railway_1000 +
  endpt_dem_tpi_500_50m + endpt_solar_radiation_10m_july +
  endpt_reservoirs_250 + endpt_lakes_250 + endpt_spraverageallyrs_500_100m_2 + I(endpt_spraverageallyrs_500_100m_2^2) +
  endpt_digestible_biomass_summer

# f <- case_ ~ strata(step_id) + sl_ + scale(startpt_roadsXXX) + sl_:scale(startpt_roadsXXX) +
#   scale(startpt_cabinsXXX) + sl_:scale(startpt_cabinsXXX) +
#   scale(endpt_roadsXXX) + scale(endpt_cabinsXXX)
f <- case_ ~ strata(step_id) + sl_ + startpt_roadsXXX + sl_:startpt_roadsXXX +
  startpt_cabinsXXX + sl_:startpt_cabinsXXX +
  endpt_roadsXXX + endpt_cabinsXXX
f <- add_zoi_formula(f, zoi_radius = c(1000, 2000, 3000), pattern = "XXX")

# data <- dat
i <- 0

sampling_params <- set_sampling_params(max_validation_size_blockH0 = 100,
                                       max_number_blockH1 = 10,
                                       max_fit_size = 1000)

sampling_params <- set_sampling_params(max_validation_size_blockH0 = 1000,
                                       max_number_blockH1 = 10,
                                       max_fit_size = 1000)

issf <- spat_strat_fit_issf(f, data, spStrat, sampling_params, i=0, kernel_vars = "sl_")
issf
hist((issf$d_validation+1)/2, xlim=c(0,1), breaks=20)
abline(v=0.5, col="black", lty="dotted", lwd=2) # <: No discrimination
abline(v=0.7, col="red", lty="solid", lwd=2) # <: poor discrimination, >: acceptable
abline(v=0.8, col="green", lty="solid", lwd=2) # >: Excellent

(issf$validation_score+1)/2
#	0.5 = No discrimination
#	0.5-0.7 = Poor discrimination
#	0.7-0.8 = Acceptable discrimination
#	0.8-0.9= Excellent discrimination
#	>0.9 = Outstanding discrimination

# in parallel ------------------------------

spat_strat_issf(f, data, spStrat, sampling_params, i=0,  kernel_vars = "sl_", out_dir="/data/R/Prosjekter/Rein/test/")

library(parallel)
library(doParallel)
library(foreach)

(starttime <- Sys.time())

cl <- makeCluster(5)
registerDoParallel(cl)
#tmp <- foreach(i=c(1:20)) %dopar% spat_strat_issf(f, data, spStrat, sampling_params, i, out_dir="/data/R/Prosjekter/Rein/test/")
tmp <- foreach(i=c(1:5)) %dopar% spat_strat_fit_issf(f, data, spStrat, sampling_params, i,
                                                     kernel_vars = "sl_",
                                                     out_dir = "/data/P-Prosjekter/41203800_oneimpact/04_tools/support_oneimpact/")
stopCluster(cl)

(endtime <- Sys.time()-starttime) # about 48.14381 mins

# Multi-model inference ------------------------------
i <- list.files("/data/R/Prosjekter/Rein/test/")
i <- i[grepl("spat_strat_issf_i", i)]
i <- gsub("spat_strat_issf_i", "", i)
i <- gsub(".rda", "", i)
i

score2weight <- function(x){
  x <- x$habitat_validation_score
  #x <- 1/mean(1/x)
  x <- min(x)
  x <- ifelse(x<0.7, 0.7, x) #truncate poorly validated models
  return(x)
}

missf <- multi_issf(i, f, data, score2weight=score2weight, weights_function=NULL, out_dir="/data/R/Prosjekter/Rein/test/")
save(missf, file="/data/R/Prosjekter/Rein/test/missf.rda")

str(missf)
(missf$data_summary[c(1:11), c(1:10)])

#weighted validation
apply(missf$validation_score,2,mean) %*% missf$weights #very good
apply(missf$validation_score,2,min) %*% missf$weights #very good

#weighted validation based on habitat only
apply(missf$habitat_validation_score,2,mean) %*% missf$weights #very good
apply(missf$habitat_validation_score,2,min) %*% missf$weights #very good
apply(missf$habitat_validation_score,2,quantile, probs=0.1) %*% missf$weights #very good
missf$habitat_validation_score

cbind(apply(missf$habitat_validation_score,2,min) , missf$weights)


hist(missf$weights, breaks=10)

sum(missf$weights)
