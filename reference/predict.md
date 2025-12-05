# Prediction of a bag of models to new data

The function `predict` makes a prediction for new data based wither on a
bag of models or on its formula, coefficients, and weights. The
prediction can be made either for a complete new dataset with all the
variables included in the formula or to predict the specific response on
one single or a group of variables in the model. In this case, all the
other variables are set to their median or mean value, to to zero
(defined by the `baseline` parameter). What controls that is which
columns are added in the `newdata` data.frame.

## Usage

``` r
predict(x, newdata, ...)

# S3 method for class 'bag'
predict(
  x,
  newdata,
  data = NULL,
  type = c("linear", "exponential", "exp", "logit", "cloglog")[1],
  wmean = TRUE,
  wq_probs = NULL,
  include = "all",
  baseline = c("median", "mean", "zero")[1],
  zoi = FALSE,
  zoi_shape = c("exp_decay", "gaussian_decay", "linear_decay", "threshold_decay")[1],
  which_cumulative = "cumulative",
  type_feature = c("point", "line", "area")[1],
  type_feature_recompute = FALSE,
  n_features = 1,
  zoi_limit = 0.05,
  resolution = 100,
  line_value = 1,
  ...
)

# S3 method for class 'formula'
predict(
  x,
  newdata,
  coefs,
  weights = 1,
  type = c("linear", "exponential", "exp", "logit", "cloglog")[1],
  wmean = TRUE,
  wq_probs = NULL,
  include = "all",
  ...
)
```

## Arguments

- x:

  `[bag,list or formula]`  
  A bag of models, resulting from a call to
  [`bag_models()`](https://ninanor.github.io/oneimpact/reference/bag_models.md),
  or a `formula` used to fit the models in the bag.

- newdata:

    
  New data set to be used for prediction. It can include all the
  variables in the formula or only those for which the user is
  interested in making a prediction from.

- ...:

    
  Additional parameters. None implemented.

- type:

  `[character="linear"]{"linear", "exponential", "exp", "logit", "cloglog"}`  
  Type of prediction. One of `"linear"` (default), `"exp"` or
  `"exponential"`, `"logit"`, or `"cloglog"`.

- wmean:

  `[logical=TRUE]`  
  Should the weighted mean values be predicted? Default is `TRUE`.

- wq_probs:

  `[vector,numeric(3)=c(0.025, 0.5, 0.975)]`  
  A three element vector with lower, mid, and higher weighted quantiles
  to be computed.

- include:

  `[character="all"]`  
  String of vector of strings with the terms (or unique parts of terms)
  to be predicted for. This does not restrict which terms we are
  focusing on - this is done by the definition of the `newdata` dataset
  and by which columns are in there. What the `include` parameters does
  is to set which other variables will be used for prediction, at their
  mean or median values, for instance.

- baseline:

  `[character="median"]{"median", "mean", "zero")}`  
  What values to choose for the baseline, i.e., for all other
  variables/terms not contained in `newdata`. It can be one of `median`,
  `"mean"`, or `"zero"`.

- zoi:

  `[logical(1)=FALSE]`  
  Are the columns in `newdata` supposed to represent zones of influence
  (ZOI) variables? This parameter should be set to `TRUE` if you
  provided a set of distances from a source that need to be translated
  into ZOI variables (cumulative or nearest ZOI from sources).

- zoi_shape:

  `[character="exp_decay"]{"exp_decay", "gaussian_decay", "linear_decay", "threshold_decay"}`  
  Shape of the zone of influence (ZOI), if `zoi = TRUE`. Default is
  `exp_decay"`. It can assume any of the possible values for the
  argument `type` in the function
  [`dist_decay()`](https://ninanor.github.io/oneimpact/reference/zoi_functions.md).

- which_cumulative:

  `[character="cumulative"]`  
  Which string or pattern to be searched on the column names of
  `newdata` and on the original data used to fit the models to represent
  the cumulative ZOI. It is used to break the names of the columns/terms
  in the formula and get the ZOI radii as numbers, to be able to create
  all the ZOI radii included in the model or bag of models.

- type_feature:

  `[character="point"]{"point", "line", "area"}`  
  Type of feature we are predicting for, for zone of influence-type
  variables. Default is `"point"`. If `type_feature = "line"`, a line is
  simulated with the function
  [`create_linear_feature_zoi()`](https://ninanor.github.io/oneimpact/reference/create_linear_feature_zoi.md)
  to get the values and account for the number of pixels of each single
  linear feature in the neighborhhod and correclty estimate the effect
  of each linear feature ZOI. The option `"area"` is still not
  implemented and for now is treated as a point feature at the origin.

- n_features:

  `[numeric(1)=1]`  
  Number of features to be used for prediction, for ZOI variables.
  Default is 1.

- resolution:

  `[numeric(1)=100]`  
  Resolution for the raster created in `create_line_feature_zoi()`, when
  `type_feature = "line"`.

- line_value:

  `[numeric(1)=1]`  
  Value set to the raster line created by `create_line_feature_zoi()`,
  when `type_feature = "line"`. It could be changed to different values
  if we want to represent e.g. the value in the linear feature as the
  roads traffic or another value for spatio-temporally dynamic
  variables.

- coefs:

  `[vector,numeric]`  
  Either a named vector of coefficients (in case there is only one
  model) or a matrix of coefficients, with rownames as the term names
  and columns as the different models/resamples. Only relevant if `x` is
  a formula.

- weights:

  `[vector,numeric=1]`  
  Vector of weights for the different models/resamples, i.e. the column
  from the `coefs` object with coefficients. A single number (by
  default, 1) in case there is only one model (`coefs` is a vector).
  Only relevant if `x` is a formula.

## See also

[`plot_response()`](https://ninanor.github.io/oneimpact/reference/plot_response.md),
`create_line_feature_zoi()`.

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
#> [1] "Starting random sampling..."

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
#> Error in UseMethod("predict"): no applicable method for 'predict' applied to an object of class "formula"

# whole bag, weighted mean - here all weights = 1
predict(x = f,
        newdata = dfvar,
        coefs = bag_object$coef,
        include = names(dfvar))
#> Error in UseMethod("predict"): no applicable method for 'predict' applied to an object of class "formula"

# whole bag, for each model separately
bag_predict(x = f,
            newdata = dfvar,
            coefs = bag_object$coef,
            wmean = FALSE,
            include = names(dfvar))
#> Error in bag_predict(x = f, newdata = dfvar, coefs = bag_object$coef,     wmean = FALSE, include = names(dfvar)): could not find function "bag_predict"

#---
# prediction using bag

# prediction for the very same dataset, linear scale
predict(x = bag_object,
        newdata = dat,
        data = dat)
#> Error in UseMethod("predict"): no applicable method for 'predict' applied to an object of class "c('bag', 'list')"

# non ZOI variable
# new data, looking only at PCA3
dfvar = data.frame(norway_pca_klima_axis3 = seq(min(bag_object$data_summary$norway_pca_klima_axis3),
                                                max(bag_object$data_summary$norway_pca_klima_axis3),
                                                length.out = 100))

predict(x = bag_object,
        newdata = dfvar,
        data = dat)
#> Error in UseMethod("predict"): no applicable method for 'predict' applied to an object of class "c('bag', 'list')"

# ZOI variable
# new data, looking only at private cabins
dfvar = data.frame(private_cabins = 1e3*seq(0.2, 20, length.out = 100))

# prediction for 1 feature, linear scale
predict(x = bag_object,
        newdata = dfvar,
        data = dat,
        zoi = TRUE,
        baseline = "zero")
#> Error in UseMethod("predict"): no applicable method for 'predict' applied to an object of class "c('bag', 'list')"

# prediction for 30 features, exp scale, with weighted confidence intervals
predict(x = bag_object,
        newdata = dfvar,
        data = dat,
        type = "exp",
        wq_probs = c(0.025, 0.975),
        zoi = TRUE,
        n_features = 30,
        baseline = "zero")
#> Error in UseMethod("predict"): no applicable method for 'predict' applied to an object of class "c('bag', 'list')"

# plot
plot(dfvar[,1],
     predict(x = bag_object,
             newdata = dfvar,
             data = dat,
             type = "exp",
             zoi = TRUE,
             n_features = 30,
             baseline = "zero")[,1])
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'y' in selecting a method for function 'plot': no applicable method for 'predict' applied to an object of class "c('bag', 'list')"
```
