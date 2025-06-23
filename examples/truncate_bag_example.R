#-------
# fit a bag

#---
# fit a bag to be tested

# load packages
library(glmnet)

# load data
data("reindeer_rsf")
# rename it just for convenience
dat <- reindeer_rsf

# formula initial structure
f <- use ~ private_cabins_cumulative_XXX + public_cabins_high_cumulative_XXX +
  trails_cumulative_XXX +
  NORUTreclass +
  # poly(norway_pca_klima_axis1, 2, raw = TRUE) +
  # poly(norway_pca_klima_axis2, 2, raw = TRUE) +
  norway_pca_klima_axis1 + norway_pca_klima_axis1_sq +
  norway_pca_klima_axis2 + norway_pca_klima_axis2_sq +
  norway_pca_klima_axis3 + norway_pca_klima_axis4

# add ZOI terms to the formula
zois <- c(100, 250, 500, 1000, 2500, 5000, 10000, 20000)
ff <- add_zoi_formula(f, zoi_radius = zois, pattern = "XXX",
                      cumulative = "",
                      type = c("exp_decay"),#, "nearest_exp_decay"),
                      separator = "", predictor_table = TRUE)
f <- ff$formula
pred_table <- ff$predictor_table

# sampling - random sampling
set.seed(1234)
samples <- create_resamples(y = dat$use,
                            p = c(0.2, 0.2, 0.2),
                            times = 20,
                            colH0 = NULL)

# fit multiple models
fittedl <- bag_fit_net_logit(f,
                             data = dat,
                             samples = samples,
                             standardize = "internal", # glmnet does the standardization of covariates
                             metric = "AUC",
                             method = "AdaptiveLasso",
                             predictor_table = pred_table,
                             parallel = "mclapply",
                             mc.cores = 8) #2)

# bag models in a single object
bag_object <- bag_models(fittedl, dat, score_threshold = 0.7)

#---------
# truncate model
bag_object_trunc <- truncate_bag(bag_object,
                             data = dat,
                             measure = "cross",
                             criterion = "first_coef",
                             wmean = FALSE)

# compare validation scores
bag_object$validation_score - bag_object_trunc$validation_score

bag_object$weighted_validation_score
bag_object_trunc$weighted_validation_score

# plot curves to check weirdness

# ZOI public cabins cumulative
dfvar = data.frame(trails_cumulative = 1e3*seq(0.2, 20, length.out = 100))

# look into curve
# plot_response(bag_object,
#               dfvar = dfvar,
#               data = dat,
#               type = "linear", zoi = TRUE,
#               type_feature =  "line",
#               type_feature_recompute = TRUE,
#               resolution = 300,
#               ci = FALSE, indiv_pred = TRUE)
# original plot
plot_response(bag_object,
              dfvar = dfvar,
              data = dat,
              type = "linear", zoi = TRUE,
              ci = FALSE, indiv_pred = TRUE)
# truncated plot
plot_response(bag_object_trunc,
              dfvar = dfvar,
              data = dat,
              type = "linear", zoi = TRUE,
              # type_feature =  "line",
              # type_feature_recompute = TRUE,
              # resolution = 300,
              ci = FALSE, indiv_pred = TRUE)

# check weirdness
weirdness(bag_object_trunc, dat)
