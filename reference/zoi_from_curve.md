# Get estimates of zone of influence (ZOI) metrics from response curves

This generic function computes ZOI metrics (maximum effect size, ZOI
radius, and impact) for ZOI predictor variables based on response curves
from statistical models. The ZOI radius is estimated as the
distance/radius at which the relative selection strength decays to a
given percentage of the maximum effect size (e.g. 95% ZOI radius for the
distance at which the effect drops to 5% of the maximum). The impact
accounts for both the effect size and the ZOI radius and corresponds to
the area under (or over, if negative) the ZOI response curve. The
function computes ZOI metrics based on from weighted summary response
curves (mean or median), and based on that computes the confidence
interval bounds or individual model ZOI metrics to represent uncertainty
in the ZOI metrics. The function supports two types of input: a
`data.frame` of individual model predictions or a `bag` object
containing an ensemble of models.

## Usage

``` r
zoi_from_curve(x, ...)

# S3 method for class 'data.frame'
zoi_from_curve(
  x,
  weights,
  percentage = 0.95,
  curve = c("median", "mean")[1],
  wq_probs = c(0.025, 0.975),
  ci = TRUE,
  type = c("linear", "exp")[1],
  mean_col_name = "mean",
  median_col_name = "quantile:0.5",
  NAasZero = TRUE
)

# S3 method for class 'bag'
zoi_from_curve(
  x,
  data,
  include = "all",
  percentage = 0.95,
  curve = c("median", "mean")[1],
  type = c("linear", "exp")[1],
  return_predictions = FALSE,
  return_format = c("list", "df")[2],
  ci = TRUE,
  wq_probs = c(0.025, 0.975),
  format_long = TRUE,
  n_features = 1,
  mean_col_name = "mean",
  median_col_name = "quantile:0.5",
  NAasZero = TRUE,
  radius_max = NULL,
  baseline = "zero",
  type_feature = "line",
  type_feature_recompute = TRUE,
  resolution = 200,
  radii = c(100, 250, 500, 1000, 2500, 5000, 10000),
  zoi_shape = c("circle", "Gauss", "rectangle", "exp_decay", "bartlett", "threshold",
    "mfilter")[1],
  ...
)
```

## Arguments

- x:

  `[data.frame,bag]`  
  Either a `data.frame` containing response curve predictions for a
  single variable, or a `bag` object containing an ensemble of models.

- ...:

  `[any]`  
  Additional arguments passed to the appropriate method.

- weights:

  `[numeric]`  
  Numeric vector of model weights used to compute the weighted mean,
  weighted median, and weighted quantiles from the individual model
  prediction columns. This should match the number of individual model
  prediction columns in `x`.

- percentage:

  `[numeric(1)=0.95]`  
  Numeric between 0 and 1. Defines the threshold for ZOI radius as a
  proportion of the maximum effect size. Default is `0.95`.

- curve:

  `[character(1)=c("mean", "median")]`  
  Character vector. Which central tendency curves to use: `"median"`,
  `"mean"`, or both.

- wq_probs:

  `[numeric,vector=c(0.025, 0.975)]`  
  Numeric vector of quantiles used to compute confidence-interval for
  the ZOI metrics when `ci = TRUE`. Ignored when `ci = FALSE`.

- ci:

  `[logical(1)=TRUE]`  
  Logical. When `TRUE`, returns ZOI metrics for weighted mean and/or
  median, and the weighted confidence interval curves. When `FALSE`,
  returns ZOI metrics for weighted mean/median and all for each
  individual model prediction curve present in the input.

- type:

  `[character(1)="linear"]{"linear", "exp"}`  
  Character. Defines whether the calculation of ZOI should be based on
  the prediction on the linear scale or the response (exponential)
  scale.

- mean_col_name:

  `[character="mean"]`  
  Name of the column containing the weighted mean response curve.

- median_col_name:

  `[character="quantile:0.5"]`  
  Name of the column containing the weighted median response curve.

- NAasZero:

  `[logical(1)=TRUE]`  
  Logical. If `TRUE`, any `NA` values in the final output are replaced
  by zero.

- data:

  `[data.frame]`  
  The original dataset used for model fitting.

- include:

  `[character="all"]`  
  Character. Either `"all"` or a regex pattern to filter selected ZOI
  variables.

- return_predictions:

  `[logical=FALSE]`  
  Logical. Whether to return the prediction curves alongside ZOI
  metrics. If `TRUE`, the output is necessarily a `list` with with all
  `predictions` and the `zoi` metrics.

- return_format:

  `[character="df"]{"list", "df"}`  
  Format of the returned ZOI metrics. Either a list of data.frames (if
  `return_format = "list"`), one for each variable, or a single
  `data.frame` (default, if `return_format = "df"`).

- format_long:

  `[logical(1)=TRUE]`  
  Logical. Whether to return the ZOI metrics in long format (with a
  `zoi_metric` column) or wide format (with separate columns for each
  metric)..

- n_features:

  `[numeric=1]`  
  Number of features used in ZOI prediction. It can a single number
  (considered the same for all ZOI variables) or a vector with the same
  number of elements as ZOI variables in the model.

- radius_max:

  `[numeric=NULL]`  
  Numeric. Maximum distance/radius to use for prediction curves. If
  `NULL` (default), the maximum value present in the bag's predictor
  table is used.

- baseline:

  `[character="zero"]`  
  Character. Baseline used in
  [`predict()`](https://ninanor.github.io/oneimpact/reference/predict.md)
  (e.g., `"zero"`).

- type_feature:

  `[character="point"]`  
  Character or vector. Type of spatial feature used in
  [`predict()`](https://ninanor.github.io/oneimpact/reference/predict.md).

- type_feature_recompute:

  `[logical=FALSE]`  
  Logical. Whether to recompute spatial features within
  [`predict()`](https://ninanor.github.io/oneimpact/reference/predict.md),
  for linear features.

- resolution:

  `[numeric=200]`  
  Integer. Resolution used in the recomuptation of ZOIs for linear
  features.

- radii:

  `[vector]`  
  Numeric vector. Radii used for ZOI modeling.

- zoi_shape:

  `[character]`  
  Character. Shape of the ZOI used in the model (e.g., `"circle"`,
  `"Gauss"`, `"exp_decay"`).

## Value

For the `data.frame` method, returns a `data.frame` with columns for
each ZOI metric (`max_effect_size`, `zoi_radius`, `effect_zoi_radius`,
`impact`). When `ci = TRUE`, the rows are the weighted `mean`, `median`,
and the lower and upper CI quantiles. When `ci = FALSE`, the rows are
`mean`, `median`, and one row per individual model prediction curve
present in the input. The table can be transformed into long format if
`format_long = TRUE`.

For the `bag` method, returns either a list or a `data.frame` of ZOI
metrics for each ZOI variable in the bag. When `ci = TRUE`, the `stats`
column contains `mean`, `median`, and the CI quantile labels. When
`ci = FALSE`, the `stats` column contains one entry per individual model
curve. If `format_long = TRUE`, the output `data.frame` is in long
format, with a `zoi_metric` column indicating the type of ZOI metric
(e.g., `max_effect_size`, `zoi_radius`, `impact`) and a `metric_value`
column with the corresponding values. If `return_predictions = TRUE`,
the function returns a list with two elements: `predictions`, which is a
list of data.frames containing the prediction curves for each ZOI
variable, and `zoi`, which contains the ZOI metrics as described above.

## See also

[`predict()`](https://ninanor.github.io/oneimpact/reference/predict.md),
[`plot_response()`](https://ninanor.github.io/oneimpact/reference/plot_response.md),
[`implausibility()`](https://ninanor.github.io/oneimpact/reference/implausibility.md)

## Examples

``` r
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
f <- use ~ cabins_private_XXX + cabins_public_XXX +
  trails_XXX +
  NORUTreclass +
  norway_pca_klima_axis1 + norway_pca_klima_axis1_sq +
  norway_pca_klima_axis2 + norway_pca_klima_axis2_sq +
  norway_pca_klima_axis3 + norway_pca_klima_axis4

# add ZOI terms to the formula
zois <- c(100, 250, 500, 1000, 2500, 5000, 10000)
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
#> [1] "Starting random sampling..."

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
dfvar <- data.frame(cabins_private = 1e3*seq(0.2, 12, length.out = 100))
pred <- predict(x = bag_object,
                newdata = dfvar,
                data = dat,
                type = "linear",
                wmean = FALSE,
                zoi = TRUE,
                n_features = 1,
                baseline = "zero")
#> Error in UseMethod("predict"): no applicable method for 'predict' applied to an object of class "c('bag', 'list')"

# df with prediction
x <- cbind(dfvar, pred)
#> Error: object 'pred' not found
(tab1 <- zoi_from_curve(x, weights = bag_object$weights, curve = "median"))
#> Error: object 'x' not found

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
           xmin = tab1$zoi_radius[2], xmax = tab1$zoi_radius[3], size = 0.5, col = "red") +
  # add median and CI for maximum effect size
  annotate("pointrange", x = 0, y = tab1$max_effect_size[1],
           ymin = tab1$max_effect_size[2], ymax = tab1$max_effect_size[3], size = 0.5,
           col = "red") +
  xlim(0, 12000)
#> Error: object 'tab1' not found

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
#> Error in UseMethod("predict"): no applicable method for 'predict' applied to an object of class "c('bag', 'list')"

# df with prediction
x <- cbind(dfvar, pred)
#> Error: object 'pred' not found
(tab30 <- zoi_from_curve(x, weights = bag_object$weights, curve = "median"))
#> Error: object 'x' not found

# plot lines
p <- plot_response(bag_object,
                   dfvar = dfvar,
                   data = dat,
                   type = "linear",
                   zoi = TRUE,
                   plot_mean = FALSE,
                   n_features = 30,
                   ci = FALSE,
                   indiv_pred = TRUE,
                   logx = FALSE)
p +
  # add median and CI for ZOI
  annotate("pointrange", x = tab30$zoi_radius[1], y = tab30$effect_zoi_radius[1],
           xmin = tab30$zoi_radius[2], xmax = tab30$zoi_radius[3], size = 0.5, col = "red") +
  # add median and CI for maximum effect size
  annotate("pointrange", x = 0, y = tab30$max_effect_size[1],
           ymin = tab30$max_effect_size[2], ymax = tab30$max_effect_size[3], size = 0.5, col = "red") +
  xlim(0, 12000)
#> Error: object 'tab30' not found

# additive effect on both max_effect_size and impact
# cbind(30*impact_of_one_cabin, impact_of_30_cabins)
tibble::tibble(zoi_measure = colnames(tab1)[c(1,4)],
               one_feature = unlist(tab1[2,c(1,4)]),
               one_feature_x30 = unlist(tab1[2,c(1,4)]*30),
               thirty_features = unlist(tab30[2,c(1,4)]))
#> Error: object 'tab1' not found

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
#> Error in UseMethod("predict"): no applicable method for 'predict' applied to an object of class "c('bag', 'list')"

# df with prediction
x <- cbind(dfvar, pred)
#> Error: object 'pred' not found
(tab_exp1 <- zoi_from_curve(x, weights = bag_object$weights, type = "exp"))
#> Error: object 'x' not found

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
  xlim(0, 12000)
#> Error: object 'tab_exp1' not found

#----
# # compute ZOI for 30 features - exponential response
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
#> Error in UseMethod("predict"): no applicable method for 'predict' applied to an object of class "c('bag', 'list')"

# df with prediction
x <- cbind(dfvar, pred)
#> Error: object 'pred' not found
(tab_exp30 <- zoi_from_curve(x, weights = bag_object$weights, type = "exp"))
#> Error: object 'x' not found

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
  xlim(0, 12000)
#> Error: object 'tab_exp30' not found

# max effect size is multiplicative (power), but not impact
# cbind(30*impact_of_one_cabin, impact_of_30_cabins)
tibble::tibble(zoi_measure = colnames(tab_exp1)[c(1,4)],
               one_feature = unlist(tab_exp1[2,c(1,4)]),
               one_feature_x30 = unlist(tab_exp1[2,c(1,4)]**30),
               thirty_features = unlist(tab_exp30[2,c(1,4)]))
#> Error: object 'tab_exp1' not found

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
                zoi_vals = c(100, 250, 500, 1000, 2500, 5000, 10000),
                resolution = 200)
#> Error in UseMethod("predict"): no applicable method for 'predict' applied to an object of class "c('bag', 'list')"

# df with prediction
x <- cbind(dfvar, pred)
#> Error: object 'pred' not found
(tab_exp1_line <- zoi_from_curve(x,
                                 weights = bag_object$weights,
                                 type = "linear"))
#> Error: object 'x' not found

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
                   # zoi_vals = c(100, 250, 500, 1000, 2500, 5000, 10000),
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
  xlim(0, 12000)
#> Error: object 'tab_exp1_line' not found

#------
# Important!
# Compare summary output vs. individual-model output
(tab_ci <- zoi_from_curve(x, weights = bag_object$weights, ci = TRUE))
#> Error: object 'x' not found
(tab_models <- zoi_from_curve(x, weights = bag_object$weights, ci = FALSE))
#> Error: object 'x' not found

#------------------------------
# checking with predictions for each model, for all ZOI variables

# zoi_from curve applied to a whole bag of models

# first let's check all the fitted ZOI response curves
vars <- c("cabins_private", "cabins_public", "trails")
type_feat <- c("point", "point", "line")
rad <- unique(bag_object$parms$predictor_table$zoi_radius); rad <- rad[!is.na(rad)]
plots <- lapply(seq_along(vars), function(i) {
  df <- data.frame(col = 1e3*seq(0.002, 12.002, length.out = 1001))
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
                # zoi_vals = rad,
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
#> Error in zoi_radius[id] <- zoi_radius_main: replacement has length zero
zois
#> [1]   100   250   500  1000  2500  5000 10000

# plot
i <- 1
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
#> Error in zois$variable: $ operator is invalid for atomic vectors
print(pp + ggtitle(var))
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'print': object 'pp' not found

i <- 2
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
#> Error in zois$variable: $ operator is invalid for atomic vectors
print(pp + ggtitle(var))
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'print': object 'pp' not found

i <- 3
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
#> Error in zois$variable: $ operator is invalid for atomic vectors
print(pp + ggtitle(var))
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'print': object 'pp' not found

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
#> Error in zoi_radius[id] <- zoi_radius_main: replacement has length zero
zois
#> [1]   100   250   500  1000  2500  5000 10000

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
#> Error in UseMethod("filter"): no applicable method for 'filter' applied to an object of class "c('double', 'numeric')"
```
