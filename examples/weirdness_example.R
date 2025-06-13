#-------
# weirdness for vector of coefficients

# weirdness for coefficients for one type of ZOI variable

# set coefficients
coefs <- c(-1, -0.5, -0.1, 0.8, 0.3, -0.1)
expected_sign <- -1
weirdness(coefs, expected_sign = expected_sign)
weirdness(coefs, expected_sign = expected_sign, which_coef = "sum")

#-------
# weirdness for data.frame with (x,y) for line

# checking for lines crossing zero
x <- seq(0, 10, 0.01)
y <- -8 + 10 * x - 1.5 * x**2
df <- data.frame(x = x, y = y)
plot(x, y); abline(h = 0, col = "red")

# n crosses
weirdness(df, response = "y", measure = "n_crosses")
# auc on the opposite side of the expected sign
weirdness(df, response = "y", measure = "response_auc_opposite")
# ratio between auc above and auc on the expected sign
weirdness(df, response = "y", measure = "response_auc_ratio")

#-------
# weirdness for bag

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

# bag_object$coef %*% bag_object$weights
# sapply(fittedl, function(x) x$train_score)

#---
# plot to check

# ZOI public cabins cumulative
dfvar = data.frame(trails_cumulative = 1e3*seq(0.2, 20, length.out = 100))

# look into curve
plot_response(bag_object,
              dfvar = dfvar,
              data = dat,
              type = "linear", zoi = TRUE,
              type_feature =  "line",
              type_feature_recompute = TRUE,
              resolution = 300,
              ci = FALSE, indiv_pred = TRUE)
# with no line, just as an example
plot_response(bag_object,
              dfvar = dfvar,
              data = dat,
              type = "linear", zoi = TRUE,
              ci = FALSE, indiv_pred = TRUE)

# we try the function with the curve above, but then test how we could work with the more correct one below
# weirdness measures
weirdness(bag_object,
          data = dat)

# for each individual model
weirdness(bag_object,
          data = dat,
          wmean = FALSE)

