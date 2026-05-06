#---
# fit a bag to be tested

######## MAYBE USE A CONSTRAINED MODEL

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
ff <- add_zoi_formula(f, zoi_radius = zois, pattern = "XXX",
                      type = c("cumulative_exp_decay"),
                      separator = "", predictor_table = TRUE)
f <- ff$formula
pred_table <- ff$predictor_table

# sampling - random sampling
set.seed(1234)
samples <- create_resamples(y = dat$use,
                            p = c(0.2, 0.2, 0.2),
                            times = 10,
                            colH0 = NULL)

# set upper limit for the ZOI variables for an ecology-constrained model
# upper_limits <- ifelse(pred_table$is_zoi, 0, Inf)
upper_limits <- c(Inf, c(rep(0, 24), rep(Inf, 21))) # needs to correspond to the terms in the model matrix

# fit model
fittedl <- bag_fit_net_logit(f,
                             data = dat,
                             samples = samples,
                             standardize = "internal", # glmnet does the standardization of covariates
                             metric = "AUC",
                             method = "Lasso",
                             predictor_table = pred_table,
                             upper.limit = upper_limits,
                             parallel = "mclapply",
                             mc.cores = 2)

# bag models in a single object
bag_object <- bag_models(fittedl, dat, score_threshold = 0.7)

# bag_object <- truncate_bag(bag_object, dat)

#----
# test for private cabins

#-----
# compute ZOI using data.frames with prediction

#---
# compute ZOI for 1 feature

# prediction for each model in the bag, using wmean = FALSE

# private cabins
dfvar <- data.frame(private_cabins = 1e3*seq(0.2, 20, length.out = 100))
pred <- predict(x = bag_object,
                newdata = dfvar,
                data = dat,
                type = "linear",
                wmean = FALSE,
                zoi = TRUE,
                n_features = 1,
                baseline = "zero")

# df with prediction
x <- cbind(dfvar, pred)
(tab1 <- zoi_from_curve(x, weights = bag_object$weights, curve = "median"))

# plot lines
p <- plot_response(bag_object,
                   dfvar = dfvar,
                   data = dat,
                   type = "linear",
                   plot_mean = FALSE,
                   zoi = TRUE,
                   ci = FALSE,
                   indiv_pred = TRUE,
                   logx = FALSE)
p +
  # add median and CI for ZOI
  annotate("pointrange", x = tab1$zoi_radius[1], y = tab1$effect_zoi_radius[1],
           xmin = tab1$zoi_radius[1], xmax = tab1$zoi_radius[2], size = 0.5) +
  # add median and CI for maximum effect size
  annotate("pointrange", x = 0, y = tab1$max_effect_size[1],
           ymin = tab1$max_effect_size[1], ymax = tab1$max_effect_size[2], size = 0.5) +
  xlim(0, 5000)

#----
# # compute ZOI for 30 features
# At the linear scale, the estimated ZOI radius does not change

# prediction for each model in the bag, using wmean = FALSE

pred <- predict(x = bag_object,
                newdata = dfvar,
                data = dat,
                type = "linear",
                wmean = FALSE,
                zoi = TRUE,
                n_features = 30,
                baseline = "zero")

# df with prediction
x <- cbind(dfvar, pred)
(tab30 <- zoi_from_curve(x, weights = bag_object$weights))

# plot lines
p <- plot_response(bag_object,
                   dfvar = dfvar,
                   data = dat,
                   type = "linear",
                   zoi = TRUE,
                   n_features = 30,
                   ci = FALSE,
                   indiv_pred = TRUE,
                   logx = FALSE)
p +
  # add median and CI for ZOI
  annotate("pointrange", x = tab30$zoi_radius[1], y = tab30$effect_zoi_radius[1],
           xmin = tab30$zoi_radius[1], xmax = tab30$zoi_radius[2], size = 0.5) +
  # add median and CI for maximum effect size
  annotate("pointrange", x = 0, y = tab30$max_effect_size[1],
           ymin = tab30$max_effect_size[1], ymax = tab30$max_effect_size[2], size = 0.5) +
  xlim(0, 5000)

# additive effect on both max_effect_size and impact
# cbind(30*impact_of_one_cabin, impact_of_30_cabins)
tibble::tibble(zoi_measure = colnames(tab1)[c(1,4)],
               one_feature = unlist(tab1[2,c(1,4)]),
               one_feature_x30 = unlist(tab1[2,c(1,4)]*30),
               thirty_features = unlist(tab30[2,c(1,4)]))

#----
# compute ZOI for 1 feature - exponential response

# prediction for each model in the bag, using wmean = FALSE
pred <- predict(x = bag_object,
                newdata = dfvar,
                data = dat,
                type = "exp",
                wmean = FALSE,
                zoi = TRUE,
                n_features = 1,
                baseline = "zero")

# df with prediction
x <- cbind(dfvar, pred)
(tab_exp1 <- zoi_from_curve(x, weights = bag_object$weights, type = "exp"))

# plot lines
p <- plot_response(bag_object,
                   dfvar = dfvar,
                   data = dat,
                   type = "exp",
                   zoi = TRUE,
                   ci = FALSE,
                   indiv_pred = TRUE,
                   logx = FALSE)
p +
  # add median and CI for ZOI
  annotate("pointrange", x = tab_exp1$zoi_radius[1], y = tab_exp1$effect_zoi_radius[1],
           xmin = tab_exp1$zoi_radius[2], xmax = tab_exp1$zoi_radius[3], size = 0.5) +
  # add median and CI for maximum effect size
  annotate("pointrange", x = 0, y = tab_exp1$max_effect_size[1],
           ymin = tab_exp1$max_effect_size[2], ymax = tab_exp1$max_effect_size[3], size = 0.5) +
  xlim(0, 20000)

#----
# # compute ZOI for 30 features - expoenential response
# At the exponential scale, the estimated ZOI radius does change

# prediction for each model in the bag, using wmean = FALSE
pred <- predict(x = bag_object,
                newdata = dfvar,
                data = dat,
                type = "exp",
                wmean = FALSE,
                zoi = TRUE,
                n_features = 30,
                baseline = "zero")

# df with prediction
x <- cbind(dfvar, pred)
(tab_exp30 <- zoi_from_curve(x, weights = bag_object$weights, type = "exp"))

# plot lines
p <- plot_response(bag_object,
                   dfvar = dfvar,
                   data = dat,
                   type = "exp",
                   zoi = TRUE,
                   n_features = 30,
                   ci = FALSE,
                   indiv_pred = TRUE,
                   logx = FALSE)
p +
  # add median and CI for ZOI
  annotate("pointrange", x = tab_exp30$zoi_radius[1], y = tab_exp30$effect_zoi_radius[1],
           xmin = tab_exp30$zoi_radius[2], xmax = tab_exp30$zoi_radius[3], size = 0.5) +
  # add median and CI for maximum effect size
  annotate("pointrange", x = 0, y = tab_exp30$max_effect_size[1],
           ymin = tab_exp30$max_effect_size[2], ymax = tab_exp30$max_effect_size[3], size = 0.5) +
  xlim(0, 20000)

# max effect size is multiplicative (power), but not impact
# cbind(30*impact_of_one_cabin, impact_of_30_cabins)
tibble::tibble(zoi_measure = colnames(tab_exp1)[c(1,4)],
               one_feature = unlist(tab_exp1[2,c(1,4)]),
               one_feature_x30 = unlist(tab_exp1[2,c(1,4)]*30),
               thirty_features = unlist(tab_exp30[2,c(1,4)]))

#----
# now we use trails, a linear feature for the example
# prediction for 1 trail at exponential response scale

# prediction
dfvar = data.frame(trails_cumulative = 1e3*seq(0.2, 20, length.out = 100))
pred <- predict(x = bag_object,
                newdata = dfvar,
                data = dat,
                type = "linear",
                wmean = FALSE,
                zoi = TRUE,
                n_features = 1,
                baseline = "zero",
                type_feature = "line",
                type_feature_recompute = TRUE,
                zoi_vals = c(100, 250, 500, 1000, 2500, 5000, 10000, 20000),
                resolution = 200)

# df with prediction
x <- cbind(dfvar, pred)
(tab_exp1_line <- zoi_from_curve(x,
                                 weights = bag_object$weights,
                                 type = "linear"))

# plot
p <- plot_response(bag_object,
                   dfvar = dfvar,
                   data = dat,
                   type = "linear",
                   zoi = TRUE,
                   n_features = 1,
                   baseline = "zero",
                   type_feature = "line",
                   type_feature_recompute = TRUE,
                   zoi_vals = c(100, 250, 500, 1000, 2500, 5000, 10000, 20000),
                   resolution = 200,
                   ci = FALSE,
                   indiv_pred = TRUE,
                   logx = FALSE)
p +
  # add median and CI for ZOI
  annotate("pointrange", x = tab_exp1_line$zoi_radius[1], y = tab_exp1_line$effect_zoi_radius[1],
           xmin = tab_exp1_line$zoi_radius[2], xmax = tab_exp1_line$zoi_radius[3], size = 0.5) +
  # add median and CI for maximum effect size
  annotate("pointrange", x = 0, y = tab_exp1_line$max_effect_size[1],
           ymin = tab_exp1_line$max_effect_size[2], ymax = tab_exp1_line$max_effect_size[3], size = 0.5) +
  xlim(0, 20000)

#------
# Important!
# Compare summary output vs. individual-model output
(tab_ci <- zoi_from_curve(x, weights = bag_object$weights, ci = TRUE))
(tab_models <- zoi_from_curve(x, weights = bag_object$weights, ci = FALSE))

#------------------------------
# checking with predictions for each model, for all ZOI variables

# zoi_from curve applied to a whole bag of models

# first let's check all the fitted ZOI response curves
vars <- c("private_cabins", "public_cabins", "trails")
type_feat <- c("point", "point", "line")
rad <- unique(bag_object$parms$predictor_table$zoi_radius); rad <- rad[!is.na(rad)]
plots <- lapply(seq_along(vars), function(i) {
  df <- data.frame(col = 1e3*seq(0.002, 20.002, length.out = 1001))
  names(df) <- vars[i]
  plot_response(bag_object,
                dfvar = df,
                data = dat,
                type = "linear",
                zoi = TRUE,
                n_features = 1,
                ci = FALSE,
                indiv_pred = TRUE,
                logx = FALSE,
                type_feature = type_feat[i],
                type_feature_recompute = TRUE,
                zoi_vals = rad,
                resolution = 200)#,
})
# plots

# we can also evaluate if there are weirdness in these curves, as
# we see them
# implausibility(bag_object,
#                data = dat,
#                type_feature = c("point", "point", "line"))

#---
# now we compute ZOI parameters for all variables, with mean, median and CI
zois <- zoi_from_curve(x = bag_object,
                       data = dat,
                       type = "linear",
                       wq_probs = c(0.025, 0.975),
                       ci = TRUE,
                       n_features = 1,
                       baseline = "zero",
                       type_feature = c("point", "point", "line"),
                       type_feature_recompute = TRUE,
                       resolution = 200,
                       zoi_shape = "exp_decay")
zois

# plot
i <- 1
# for(i in seq_along(vars)) {
var <- vars[i]
pp <- plots[[i]] +
  annotate("pointrange",
           x = zois[grepl(var, zois$variable) & zois$zoi_metric == "zoi_radius" & zois$stats == "quantile:0.5",]$metric_value,
           y = zois[grepl(var, zois$variable) & zois$zoi_metric == "effect_zoi_radius" & zois$stats == "quantile:0.5",]$metric_value,
           xmin = zois[grepl(var, zois$variable) & zois$zoi_metric == "zoi_radius" & zois$stats == "quantile:0.025",]$metric_value,
           xmax = zois[grepl(var, zois$variable) & zois$zoi_metric == "zoi_radius" & zois$stats == "quantile:0.975",]$metric_value,
           size = 0.5) +
  annotate("pointrange",
           x = 0,
           y = zois[grepl(var, zois$variable) & zois$zoi_metric == "max_effect_size" & zois$stats == "quantile:0.5",]$metric_value,
           ymin = zois[grepl(var, zois$variable) & zois$zoi_metric == "max_effect_size" & zois$stats == "quantile:0.025",]$metric_value,
           ymax = zois[grepl(var, zois$variable) & zois$zoi_metric == "max_effect_size" & zois$stats == "quantile:0.975",]$metric_value,
           size = 0.5) +
  xlim(0, 5000)
print(pp + ggtitle(var))

i <- 2
# for(i in seq_along(vars)) {
var <- vars[i]
pp <- plots[[i]] +
  annotate("pointrange",
           x = zois[grepl(var, zois$variable) & zois$zoi_metric == "zoi_radius" & zois$stats == "quantile:0.5",]$metric_value,
           y = zois[grepl(var, zois$variable) & zois$zoi_metric == "effect_zoi_radius" & zois$stats == "quantile:0.5",]$metric_value,
           xmin = zois[grepl(var, zois$variable) & zois$zoi_metric == "zoi_radius" & zois$stats == "quantile:0.025",]$metric_value,
           xmax = zois[grepl(var, zois$variable) & zois$zoi_metric == "zoi_radius" & zois$stats == "quantile:0.975",]$metric_value,
           size = 0.5) +
  annotate("pointrange",
           x = 0,
           y = zois[grepl(var, zois$variable) & zois$zoi_metric == "max_effect_size" & zois$stats == "quantile:0.5",]$metric_value,
           ymin = zois[grepl(var, zois$variable) & zois$zoi_metric == "max_effect_size" & zois$stats == "quantile:0.025",]$metric_value,
           ymax = zois[grepl(var, zois$variable) & zois$zoi_metric == "max_effect_size" & zois$stats == "quantile:0.975",]$metric_value,
           size = 0.5)
print(pp + ggtitle(var))

i <- 3
# for(i in seq_along(vars)) {
var <- vars[i]
pp <- plots[[i]] +
  annotate("pointrange",
           x = zois[grepl(var, zois$variable) & zois$zoi_metric == "zoi_radius" & zois$stats == "quantile:0.5",]$metric_value,
           y = zois[grepl(var, zois$variable) & zois$zoi_metric == "effect_zoi_radius" & zois$stats == "quantile:0.5",]$metric_value,
           xmin = zois[grepl(var, zois$variable) & zois$zoi_metric == "zoi_radius" & zois$stats == "quantile:0.025",]$metric_value,
           xmax = zois[grepl(var, zois$variable) & zois$zoi_metric == "zoi_radius" & zois$stats == "quantile:0.975",]$metric_value,
           size = 0.5) +
  annotate("pointrange",
           x = 0,
           y = zois[grepl(var, zois$variable) & zois$zoi_metric == "max_effect_size" & zois$stats == "quantile:0.5",]$metric_value,
           ymin = zois[grepl(var, zois$variable) & zois$zoi_metric == "max_effect_size" & zois$stats == "quantile:0.025",]$metric_value,
           ymax = zois[grepl(var, zois$variable) & zois$zoi_metric == "max_effect_size" & zois$stats == "quantile:0.975",]$metric_value,
           size = 0.5)
print(pp + ggtitle(var))

#-----------
# Now we can also compute the values for each model and look at the distribution of ZOI values,
# beyond the confidence interval

# computing zoi metrics for each model, using ci = FALSE
zois <- zoi_from_curve(x = bag_object,
                       data = dat,
                       type = "linear",
                       ci = FALSE,
                       n_features = 1,
                       baseline = "zero",
                       type_feature = c("point", "point", "line"),
                       type_feature_recompute = TRUE,
                       resolution = 200,
                       zoi_shape = "exp_decay")
zois

# making a violin plot based on the distribution of ZOI values
zois |>
  dplyr::filter(grepl("Resample", stats), zoi_metric == "zoi_radius") |>
  ggplot(aes(x = variable, y = metric_value)) +
  geom_violin(fill = "red", alpha = 0.3) +
  coord_flip() +
  labs(y = "ZOI radius",
       x = "Disturbance") +
  # scale_y_log10() +
  theme_minimal()
