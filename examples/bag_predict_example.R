#---
# fit a bag to be tested

# load packages
library(glmnet)
library(ggplot2)

# load data
data("reindeer_rsf")
# rename it just for convenience
dat <- reindeer_rsf

# formula initial structure
f <- use ~ private_cabins_XXX + public_cabins_high_XXX +
  trails_XXX +
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
                     separator = "", predictor_table = TRUE)$formula

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
# prediction using formula

# new data, looking only at PCA1
dfvar = data.frame(norway_pca_klima_axis1 = seq(min(bag_object$data_summary$norway_pca_klima_axis1),
                                                max(bag_object$data_summary$norway_pca_klima_axis1),
                                                length.out = 100))
dfvar$norway_pca_klima_axis1_sq = dfvar$norway_pca_klima_axis1**2

# one model only
predict(x = f,
        newdata = dfvar,
        coefs = bag_object$coef[,1],
        include = "axis1")

# whole bag, weighted mean - here all weights = 1
predict(x = f,
        newdata = dfvar,
        coefs = bag_object$coef,
        include = names(dfvar))

# whole bag, for each model separately
bag_predict(x = f,
            newdata = dfvar,
            coefs = bag_object$coef,
            wmean = FALSE,
            include = names(dfvar))

#---
# prediction using bag

# prediction for the very same dataset, linear scale
predict(x = bag_object,
        newdata = dat,
        data = dat)

# non ZOI variable
# new data, looking only at PCA3
dfvar = data.frame(norway_pca_klima_axis3 = seq(min(bag_object$data_summary$norway_pca_klima_axis3),
                                                max(bag_object$data_summary$norway_pca_klima_axis3),
                                                length.out = 100))

predict(x = bag_object,
        newdata = dfvar,
        data = dat)

# ZOI variable
# new data, looking only at private cabins
dfvar = data.frame(private_cabins = 1e3*seq(0.2, 20, length.out = 100))

# prediction for 1 feature, linear scale
predict(x = bag_object,
        newdata = dfvar,
        data = dat,
        zoi = TRUE,
        baseline = "zero")

# prediction for 30 features, exp scale, with weighted confidence intervals
predict(x = bag_object,
        newdata = dfvar,
        data = dat,
        type = "exp",
        wq_probs = c(0.025, 0.975),
        zoi = TRUE,
        n_features = 30,
        baseline = "zero")

# plot
plot(dfvar[,1],
     predict(x = bag_object,
             newdata = dfvar,
             data = dat,
             type = "exp",
             zoi = TRUE,
             n_features = 30,
             baseline = "zero")[,1])
