# Predict from a bag of models or from formula + coefficients

This function predicts responses for new data using either:

- a `bag` object (via `predict.bag`), or

- a model `formula` plus `coefs` and `weights` (via `predict.formula`).

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
  A bag of models from
  [`bag_models()`](https://ninanor.github.io/oneimpact/reference/bag_models.md),
  or a `formula` used to build the model matrix for prediction.

- newdata:

  `[data.frame]`  
  New data for prediction. Can contain all model variables or only focal
  variables.

- ...:

  `[any]`  
  Additional arguments passed to methods.

- data:

  `[data.frame=NULL]`  
  Original data used for model fitting. Used in `predict.bag()` to
  recover baseline/categorical levels. Not used by `predict.formula()`.

- type:

  `[character="linear"]{"linear", "exponential", "exp", "logit", "cloglog"}`  
  Prediction scale.

- wmean:

  `[logical(1)=TRUE]`  
  Whether to compute and return weighted mean prediction.

- wq_probs:

  `[numeric,vector=c(0.025, 0.5, 0.975)]`  
  Weighted quantile probabilities for prediction summaries. If `NULL`,
  quantiles are not returned.

- include:

  `[character="all"]`  
  Terms to include in prediction. Use `"all"` or one/more string
  patterns matching term names.

- baseline:

  `[character="median"]{"median", "mean", "zero"}`  
  Baseline values for non-focal predictors.

- zoi:

  `[logical(1)=FALSE]`  
  If `TRUE`, columns in `newdata` are treated as distance inputs and
  transformed into ZOI predictors.

- zoi_shape:

  `[character="exp_decay"]{"exp_decay", "gaussian_decay", "linear_decay", "threshold_decay"}`  
  ZOI decay shape used when `zoi = TRUE`.

- which_cumulative:

  `[character(1)="cumulative"]`  
  Pattern used to identify cumulative ZOI terms.

- type_feature:

  `[character="point"]{"point", "line", "area"}`  
  Feature type for ZOI prediction.

- type_feature_recompute:

  `[logical(1)=FALSE]`  
  Whether to recompute line- and area-feature geometry approximation for
  ZOI calculations.

- n_features:

  `[numeric(1)=1]`  
  Number of features used in ZOI prediction.

- zoi_limit:

  `[numeric(1)=0.05]`  
  Lower influence threshold used by non-vanishing ZOI decay functions
  (relevant for line-feature ZOI transformation). See
  [`zoi_functions()`](https://ninanor.github.io/oneimpact/reference/zoi_functions.md).

- resolution:

  `[numeric(1)=100]`  
  Raster resolution used for line-feature ZOI approximation. Used when
  recomputing the ZOI variables for line and area features.

- line_value:

  `[numeric(1)=1]`  
  Value assigned to rasterized line cells when `type_feature = "line"`.
  Used when recomputing the ZOI variables for line and area features.

- coefs:

  `[numeric vector or matrix]`  
  Coefficients used by `predict.formula()`. If matrix, rows are term
  names and columns are model/resample coefficients.

- weights:

  `[numeric=1]`  
  Model weights used by `predict.formula()` when combining predictions
  across models.

## Value

A `data.frame` (or matrix-like object) with predicted values. Output
columns depend on `wmean` and `wq_probs`:

- weighted quantiles (if `wq_probs` is not `NULL`)

- weighted mean (if `wmean = TRUE`)

- individual model predictions (if `wmean = FALSE` and
  `wq_probs = NULL`)

## Details

Predictions can be computed for full covariate data or for focal
predictors only, while non-focal predictors are set to a baseline values
(`median`, `mean`, or `zero`). It also supports ZOI-distance inputs that
are internally transformed into ZOI predictors.

## See also

[`plot_response()`](https://ninanor.github.io/oneimpact/reference/plot_response.md),
[`create_linear_feature_zoi()`](https://ninanor.github.io/oneimpact/reference/create_linear_feature_zoi.md).

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
