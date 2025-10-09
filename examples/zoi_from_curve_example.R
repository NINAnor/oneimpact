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

# fit multiple models
fittedl <- bag_fit_net_logit(f,
                             data = dat,
                             samples = samples,
                             standardize = "internal", # glmnet does the standardization of covariates
                             metric = "AUC",
                             method = "Lasso",
                             predictor_table = pred_table,
                             parallel = "mclapply",
                             mc.cores = 2)

# bag models in a single object
bag_object <- bag_models(fittedl, dat, score_threshold = 0.7)
#
# bag_object <- truncate_bag(bag_object, dat)

#----
# first we take a look at a response curve, just have an idea

# private cabins
dfvar <- data.frame(private_cabins = 1e3*seq(0.2, 20, length.out = 100))
plot_response(bag_object,
              dfvar = dfvar,
              data = dat,
              type = "exp",
              zoi = TRUE,
              ci = FALSE,
              indiv_pred = TRUE,
              logx = TRUE,
              ylim = ylim(0, 2))

#-----
# compute ZOI using data.frames with prediction

#---
# compute ZOI for 1 feature

# prediction
pred <- predict(x = bag_object,
                newdata = dfvar,
                data = dat,
                type = "linear",
                wq_probs = c(0.025, 0.5, 0.975),
                zoi = TRUE,
                n_features = 1,
                baseline = "zero")

# df with prediction
x <- cbind(dfvar, pred)
(tab1 <- zoi_from_curve(x))

# plot

# plot lines
p <- plot_response(bag_object,
                   dfvar = dfvar,
                   data = dat,
                   type = "linear",
                   zoi = TRUE,
                   ci = FALSE,
                   indiv_pred = TRUE,
                   logx = FALSE)#,
# ylim = ylim(0, 2))
p +
  annotate("pointrange", x = tab1[2,1], y = tab1[3,1],
           xmin = tab1[2,3], xmax = tab1[2,4], size = 0.5) +
  annotate("pointrange", x = 0, y = tab1[1,1],
           ymin = tab1[1,3], ymax = tab1[1,4], size = 0.5) +
  xlim(0, 5000)

#----
# # compute ZOI for 30 features
# At the linear scale, the estimated ZOI radius dos not change

# prediction
pred <- predict(x = bag_object,
                newdata = dfvar,
                data = dat,
                type = "linear",
                wq_probs = c(0.025, 0.5, 0.975),
                zoi = TRUE,
                n_features = 30,
                baseline = "zero")

# df with prediction
x <- cbind(dfvar, pred)
(tab30 <- zoi_from_curve(x))

# plot
p <- plot_response(bag_object,
                   dfvar = dfvar,
                   data = dat,
                   type = "linear",
                   zoi = TRUE,
                   n_features = 30,
                   ci = FALSE,
                   indiv_pred = TRUE,
                   logx = FALSE)#,
# ylim = ylim(0, 2))
p +
  annotate("pointrange", x = tab30[2,1], y = tab30[3,1],
           xmin = tab30[2,3], xmax = tab30[2,4], size = 0.5) +
  annotate("pointrange", x = 0, y = tab30[1,1],
           ymin = tab30[1,3], ymax = tab30[1,4], size = 0.5) +
  xlim(0, 5000)

# additive effect on both max_effect_size and impact
# cbind(30*impact_of_one_cabin, impact_of_30_cabins)
cbind(tab1[c(1,4),1]*30, tab30[c(1,4),1])

#----
# compute ZOI for 1 feature - exponential

# prediction
pred <- predict(x = bag_object,
                newdata = dfvar,
                data = dat,
                type = "exp",
                wq_probs = c(0.025, 0.5, 0.975),
                zoi = TRUE,
                n_features = 1,
                baseline = "zero")

# df with prediction
x <- cbind(dfvar, pred)
(tab_exp1 <- zoi_from_curve(x, type = "exp"))

# plot
p <- plot_response(bag_object,
                   dfvar = dfvar,
                   data = dat,
                   type = "exp",
                   zoi = TRUE,
                   n_features = 1,
                   ci = FALSE,
                   indiv_pred = TRUE,
                   logx = FALSE)#,
# ylim = ylim(0, 2))
p +
  annotate("pointrange", x = tab_exp1[2,1], y = tab_exp1[3,1],
           xmin = tab_exp1[2,3], xmax = tab_exp1[2,4], size = 0.5) +
  annotate("pointrange", x = 0, y = tab_exp1[1,1],
           ymin = tab_exp1[1,3], ymax = tab_exp1[1,4], size = 0.5) +
  xlim(0, 10000)

#----
# compute ZOI for 30 features - exponential

# prediction
pred <- predict(x = bag_object,
                newdata = dfvar,
                data = dat,
                type = "exp",
                wq_probs = c(0.025, 0.5, 0.975),
                zoi = TRUE,
                n_features = 30,
                baseline = "zero")

# df with prediction
x <- cbind(dfvar, pred)
(tab_exp30 <- zoi_from_curve(x, type = "exp"))

# plot
p <- plot_response(bag_object,
                   dfvar = dfvar,
                   data = dat,
                   type = "exp",
                   zoi = TRUE,
                   n_features = 30,
                   ci = FALSE,
                   indiv_pred = TRUE,
                   logx = FALSE)#,
# ylim = ylim(0, 2))
p +
  annotate("pointrange", x = tab_exp30[2,1], y = tab_exp30[3,1],
           xmin = tab_exp30[2,3], xmax = tab_exp30[2,4], size = 0.5) +
  annotate("pointrange", x = 0, y = tab_exp30[1,1],
           ymin = tab_exp30[1,3], ymax = tab_exp30[1,4], size = 0.5) +
  xlim(0, 10000)

# max effect size is multiplicative (power),
cbind(c(tab_exp1[c(1),1]**30, tab_exp1[c(4),1]*30),
      c(tab_exp30[c(1),1], tab_exp30[c(4),1]))

#----
# compute ZOI for 1 linear feature - exponential response
# (there are issues when it is linear, check later)

# prediction
dfvar = data.frame(trails_cumulative = 1e3*seq(0.2, 20, length.out = 100))
pred <- predict(x = bag_object,
                newdata = dfvar,
                data = dat,
                type = "exp",
                wq_probs = c(0.025, 0.5, 0.975),
                zoi = TRUE,
                n_features = 1,
                baseline = "zero",
                type_feature = "line",
                type_feature_recompute = TRUE,
                zoi_vals = c(100, 250, 500, 1000, 2500, 5000, 10000, 20000),
                resolution = 200)

# df with prediction
x <- cbind(dfvar, pred)
(tab_exp1_line <- zoi_from_curve(x, type = "exp"))

# plot
p <- plot_response(bag_object,
                   dfvar = dfvar,
                   data = dat,
                   type = "exp",
                   zoi = TRUE,
                   n_features = 1,
                   baseline = "zero",
                   type_feature = "line",
                   type_feature_recompute = TRUE,
                   zoi_vals = c(100, 250, 500, 1000, 2500, 5000, 10000, 20000),
                   resolution = 200,
                   ci = FALSE,
                   indiv_pred = TRUE,
                   logx = FALSE)#,
# ylim = ylim(0, 2))
p +
  annotate("pointrange", x = tab_exp1_line[2,1], y = tab_exp1_line[3,1],
           xmin = tab_exp1_line[2,3], xmax = tab_exp1_line[2,4], size = 0.5) +
  annotate("pointrange", x = 0, y = tab_exp1_line[1,1],
           ymin = tab_exp1_line[1,3], ymax = tab_exp1_line[1,4], size = 0.5) +
  xlim(0, 20000)

#-----
# check using linear response, then we have issues
# prediction
dfvar = data.frame(trails_cumulative = 1e3*seq(0.2, 20, length.out = 100))
pred <- predict(x = bag_object,
                newdata = dfvar,
                data = dat,
                type = "linear",
                wq_probs = c(0.025, 0.5, 0.975),
                zoi = TRUE,
                n_features = 1,
                baseline = "zero",
                type_feature = "line",
                type_feature_recompute = TRUE,
                zoi_vals = c(100, 250, 500, 1000, 2500, 5000, 10000, 20000),
                resolution = 200)

# df with prediction
x <- cbind(dfvar, pred)
(tab_exp1_line <- zoi_from_curve(x, type = "linear"))

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
                   logx = FALSE)#,
# ylim = ylim(0, 2))
p +
  annotate("pointrange", x = tab_exp1_line[2,1], y = tab_exp1_line[3,1],
           xmin = tab_exp1_line[2,3], xmax = tab_exp1_line[2,4], size = 0.5) +
  annotate("pointrange", x = 0, y = tab_exp1_line[1,1],
           ymin = tab_exp1_line[1,3], ymax = tab_exp1_line[1,4], size = 0.5) +
  xlim(0, 20000)

#------------
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
plots

# we can also evaluate if there are weirdness in these curves, as
# we see them
weirdness(bag_object,
          data = dat,
          type_feature = c("point", "point", "line"))

#---
# now we compute ZOI parameters for all variables
zois <- zoi_from_curve(x = bag_object,
                                 data = dat,
                                 type = "linear",
                                 wq_probs = c(0.025, 0.5, 0.975),
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
           x = zois[grepl(var, zois$variable) & zois$zoi_measure == "zoi_radius",]$mean,
           y = zois[grepl(var, zois$variable) & zois$zoi_measure == "effect_zoi_radius",]$mean,
           xmin = zois[grepl(var, zois$variable) & zois$zoi_measure == "zoi_radius",]$`quantile:0.025`,
           xmax = zois[grepl(var, zois$variable) & zois$zoi_measure == "zoi_radius",]$`quantile:0.975`,
           size = 0.5) +
  annotate("pointrange",
           x = 0,
           y = zois[grepl(var, zois$variable) & zois$zoi_measure == "max_effect_size",]$mean,
           ymin = zois[grepl(var, zois$variable) & zois$zoi_measure == "max_effect_size",]$`quantile:0.025`,
           ymax = zois[grepl(var, zois$variable) & zois$zoi_measure == "max_effect_size",]$`quantile:0.975`,
           size = 0.5) +
  xlim(0, 5000)
print(pp + ggtitle(var))

i <- 2
# for(i in seq_along(vars)) {
var <- vars[i]
pp <- plots[[i]] +
  annotate("pointrange",
           x = zois[grepl(var, zois$variable) & zois$zoi_measure == "zoi_radius",]$mean,
           y = zois[grepl(var, zois$variable) & zois$zoi_measure == "effect_zoi_radius",]$mean,
           xmin = zois[grepl(var, zois$variable) & zois$zoi_measure == "zoi_radius",]$`quantile:0.025`,
           xmax = zois[grepl(var, zois$variable) & zois$zoi_measure == "zoi_radius",]$`quantile:0.975`,
           size = 0.5) +
  annotate("pointrange",
           x = 0,
           y = zois[grepl(var, zois$variable) & zois$zoi_measure == "max_effect_size",]$mean,
           ymin = zois[grepl(var, zois$variable) & zois$zoi_measure == "max_effect_size",]$`quantile:0.025`,
           ymax = zois[grepl(var, zois$variable) & zois$zoi_measure == "max_effect_size",]$`quantile:0.975`,
           size = 0.5) +
  xlim(0, 5000)
print(pp + ggtitle(var))

i <- 3
# for(i in seq_along(vars)) {
var <- vars[i]
pp <- plots[[i]] +
  annotate("pointrange",
           x = zois[grepl(var, zois$variable) & zois$zoi_measure == "zoi_radius",]$mean,
           y = zois[grepl(var, zois$variable) & zois$zoi_measure == "effect_zoi_radius",]$mean,
           xmin = zois[grepl(var, zois$variable) & zois$zoi_measure == "zoi_radius",]$`quantile:0.025`,
           xmax = zois[grepl(var, zois$variable) & zois$zoi_measure == "zoi_radius",]$`quantile:0.975`,
           size = 0.5) +
  annotate("pointrange",
           x = 0,
           y = zois[grepl(var, zois$variable) & zois$zoi_measure == "max_effect_size",]$mean,
           ymin = zois[grepl(var, zois$variable) & zois$zoi_measure == "max_effect_size",]$`quantile:0.025`,
           ymax = zois[grepl(var, zois$variable) & zois$zoi_measure == "max_effect_size",]$`quantile:0.975`,
           size = 0.5) +
  xlim(0, 20000)
print(pp + ggtitle(var))
