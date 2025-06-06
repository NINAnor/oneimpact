#---
# fit a bag to be tested

# load packages
library(glmnet)

# load data
data("reindeer_rsf")
# rename it just for convenience
dat <- reindeer_rsf

# formula initial structure
f <- use ~ private_cabins_XXX + public_cabins_high_XXX +
  NORUTreclass +
  # poly(norway_pca_klima_axis1, 2, raw = TRUE) +
  # poly(norway_pca_klima_axis2, 2, raw = TRUE) +
  norway_pca_klima_axis1 + norway_pca_klima_axis1_sq +
  norway_pca_klima_axis2 + norway_pca_klima_axis2_sq +
  norway_pca_klima_axis3 + norway_pca_klima_axis4

# add ZOI terms to the formula
zois <- c(100, 250, 500, 1000, 2500, 5000, 10000, 20000)
f <- add_zoi_formula(f, zoi_radius = zois, pattern = "XXX",
                     type = c("cumulative_exp_decay"),
                     separator = "_", predictor_table = TRUE)$formula

# sampling - random sampling
set.seed(1234)
samples <- create_resamples(y = dat$use,
                            p = c(0.2, 0.2, 0.2),
                            times = 10,
                            colH0 = NULL)

# fit multiple models
fittedl <- bag_fit_net_logit(f,
                             data = dat,
                             samples = samples,
                             standardize = "internal", # glmnet does the standardization of covariates
                             metric = "AUC",
                             method = "AdaptiveLasso",
                             parallel = "mclapply",
                             mc.cores = 2)

# bag models in a single object
bag_object <- bag_models(fittedl, dat, score_threshold = 0.7)

#---
# plot predictions for non-ZOI variables

# plot for PCA1
dfvar <- data.frame(norway_pca_klima_axis1 = seq(min(bag_object$data_summary$norway_pca_klima_axis1),
                                                 max(bag_object$data_summary$norway_pca_klima_axis1),
                                                 length.out = 100))
dfvar$norway_pca_klima_axis1_sq = dfvar$norway_pca_klima_axis1**2

# plot mean response in linear scale with weighted interquartile range
plot_response(bag_object,
              dfvar = dfvar,
              data = dat,
              wq_probs = c(0.25, 0.5, 0.75),
              plot_median = FALSE) # remove median, plot only weighted mean

# plot median response in exponential scale with weighted interquartile range
plot_response(bag_object,
              dfvar = dfvar,
              data = dat,
              type = "exp",
              wq_probs = c(0.25, 0.5, 0.75),
              plot_mean = FALSE) # remove mean, plot only weighted median

#---
# plot predictions for ZOI variables

# plot for private cabins

# define newdata based only on the distances from the source (public cabins)
dfvar = data.frame(private_cabins = 1e3*seq(0.2, 20, length.out = 100))

# plot mean response in exponential scale, with individual lines, x in log scale
# in exponential scale, relative selection strength = 1 corresponse to no effect
# prediction for for 1 cabin only
plot_response(bag_object,
              dfvar = dfvar,
              data = dat,
              type = "exp",
              zoi = TRUE,
              ci = FALSE,
              indiv_pred = TRUE,
              logx = TRUE,
              ylim = ylim(0, 2))

# prediction for for 10 cabins located at the origin
plot_response(bag_object,
              dfvar = dfvar,
              data = dat,
              type = "exp",
              zoi = TRUE,
              n_features = 10,
              ci = FALSE,
              indiv_pred = TRUE,
              logx = TRUE,
              ylim = ylim(0, 2))
