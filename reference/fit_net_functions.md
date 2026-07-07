# Fits a conditional logistic regression/SSF/iSSF with penalized regression using glmnet in a cross-validation setup

This function runs one single realization of a model fit for a specific
sample composed by a fit set, a test set, and a validation set. This
function does the actual data set up and model fitting by calling
[`net_clogit()`](https://ninanor.github.io/oneimpact/reference/net_clogit.md)
and
[`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html).
The function is also used within
[`bag_fit_net_clogit()`](https://ninanor.github.io/oneimpact/reference/bag_fit_net_functions.md).

This function runs one single realization of a model fit for a specific
sample composed by a fit set, a test set, and a validation set. This
function does the actual data set up and model fitting by calling
[`net_logit()`](https://ninanor.github.io/oneimpact/reference/net_logit.md)
and
[`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html).
The function is also used within
[`bag_fit_net_logit()`](https://ninanor.github.io/oneimpact/reference/bag_fit_net_functions.md).

## Usage

``` r
fit_net_clogit(
  f,
  data,
  samples,
  i = 1,
  kernel_vars = c("step_length", "ta"),
  metric = c("coxnet.deviance", "Cindex", "conditionalAUC", "conditionalSomersD")[1],
  metrics_evaluate = c("coxnet.deviance", "Cindex", "conditionalAUC"),
  method = c("Lasso", "Ridge", "AdaptiveLasso", "DistanceDecayLasso", "DDLasso",
    "OneZOI-AdaptiveLasso", "OZ-AdaptiveLasso", "Grouped-AdaptiveLasso",
    "G-AdaptiveLasso", "ElasticNet")[1],
  alpha = NULL,
  penalty.factor = NULL,
  gamma = 1,
  standardize = c("internal", "external", FALSE)[1],
  predictor_table = NULL,
  function_lasso_decay = c(log, function(x) x/1000)[[1]],
  value_lasso_decay = 1,
  function_hypothesis = c(exp)[[1]],
  expected_sign_hypothesis = -1,
  factor_grouped_lasso = 1,
  replace_missing_NA = TRUE,
  na.action = "na.pass",
  out_dir_file = NULL,
  verbose = FALSE,
  ...
)

fit_net_ssf(
  f,
  data,
  samples,
  i = 1,
  kernel_vars = c("step_length", "ta"),
  metric = c("coxnet.deviance", "Cindex", "conditionalAUC", "conditionalSomersD")[1],
  metrics_evaluate = c("coxnet.deviance", "Cindex", "conditionalAUC"),
  method = c("Lasso", "Ridge", "AdaptiveLasso", "DistanceDecayLasso", "DDLasso",
    "OneZOI-AdaptiveLasso", "OZ-AdaptiveLasso", "Grouped-AdaptiveLasso",
    "G-AdaptiveLasso", "ElasticNet")[1],
  alpha = NULL,
  penalty.factor = NULL,
  gamma = 1,
  standardize = c("internal", "external", FALSE)[1],
  predictor_table = NULL,
  function_lasso_decay = c(log, function(x) x/1000)[[1]],
  value_lasso_decay = 1,
  function_hypothesis = c(exp)[[1]],
  expected_sign_hypothesis = -1,
  factor_grouped_lasso = 1,
  replace_missing_NA = TRUE,
  na.action = "na.pass",
  out_dir_file = NULL,
  verbose = FALSE,
  ...
)

fit_net_issf(
  f,
  data,
  samples,
  i = 1,
  kernel_vars = c("step_length", "ta"),
  metric = c("coxnet.deviance", "Cindex", "conditionalAUC", "conditionalSomersD")[1],
  metrics_evaluate = c("coxnet.deviance", "Cindex", "conditionalAUC"),
  method = c("Lasso", "Ridge", "AdaptiveLasso", "DistanceDecayLasso", "DDLasso",
    "OneZOI-AdaptiveLasso", "OZ-AdaptiveLasso", "Grouped-AdaptiveLasso",
    "G-AdaptiveLasso", "ElasticNet")[1],
  alpha = NULL,
  penalty.factor = NULL,
  gamma = 1,
  standardize = c("internal", "external", FALSE)[1],
  predictor_table = NULL,
  function_lasso_decay = c(log, function(x) x/1000)[[1]],
  value_lasso_decay = 1,
  function_hypothesis = c(exp)[[1]],
  expected_sign_hypothesis = -1,
  factor_grouped_lasso = 1,
  replace_missing_NA = TRUE,
  na.action = "na.pass",
  out_dir_file = NULL,
  verbose = FALSE,
  ...
)

fit_net_logit(
  f,
  data,
  samples,
  i = 1,
  metric = c("AUC")[1],
  metrics_evaluate = c("AUC"),
  method = c("Lasso", "Ridge", "AdaptiveLasso", "DistanceDecayLasso", "DDLasso",
    "TruncatedLasso", "OneZOI-AdaptiveLasso", "OZ-AdaptiveLasso",
    "Grouped-AdaptiveLasso", "G-AdaptiveLasso", "ElasticNet")[1],
  alpha = NULL,
  penalty.factor = NULL,
  gamma = 1,
  standardize = c("internal", "external", FALSE)[1],
  predictor_table = NULL,
  function_lasso_decay = c(log, function(x) x/1000)[[1]],
  value_lasso_decay = 1,
  function_hypothesis = c(exp)[[1]],
  expected_sign_hypothesis = -1,
  factor_grouped_lasso = 1,
  replace_missing_NA = TRUE,
  na.action = "na.pass",
  out_dir_file = NULL,
  verbose = FALSE,
  ...
)

fit_net_rsf(
  f,
  data,
  samples,
  i = 1,
  metric = c("AUC")[1],
  metrics_evaluate = c("AUC"),
  method = c("Lasso", "Ridge", "AdaptiveLasso", "DistanceDecayLasso", "DDLasso",
    "TruncatedLasso", "OneZOI-AdaptiveLasso", "OZ-AdaptiveLasso",
    "Grouped-AdaptiveLasso", "G-AdaptiveLasso", "ElasticNet")[1],
  alpha = NULL,
  penalty.factor = NULL,
  gamma = 1,
  standardize = c("internal", "external", FALSE)[1],
  predictor_table = NULL,
  function_lasso_decay = c(log, function(x) x/1000)[[1]],
  value_lasso_decay = 1,
  function_hypothesis = c(exp)[[1]],
  expected_sign_hypothesis = -1,
  factor_grouped_lasso = 1,
  replace_missing_NA = TRUE,
  na.action = "na.pass",
  out_dir_file = NULL,
  verbose = FALSE,
  ...
)

grouped_func(coefs, phi_group = 0)
```

## Arguments

- f:

  `[formula]`  
  Formula of the model to be fitted, with all possible candidate terms.

- data:

  `[data.frame,tibble]`  
  Complete data set to be analyzed.

- samples:

  `[list]`  
  List of samples with at least three elements: train, test, and
  validate. Each elements might have several elements, each representing
  the lines of `data` to be sampled for each resample. Typically, this
  is computed by the function
  [`create_resamples()`](https://ninanor.github.io/oneimpact/reference/create_resamples.md).

- i:

  `[numeric(1)=1]`  
  Index of the current resample iteration.

- kernel_vars:

  `[vector,character=c("step_length", "ta")]`  
  Vector of strings with the names of the variables related to the
  movement kernel, included in the model (for instance, `"step_length"`
  and `"turning_angle"`)

- metric:

  `[function,character]{AUC, conditionalBoyce, conditionalSomersD, conditionalAUC}`  
  Function representing the metric to evaluate goodness-of-fit. One of
  AUC (Default), conditionalBoyce, conditionalSomersD, and
  conditionalAUC. A user-defined function might be provided, with a
  condition that it must be maximized to find the best fit model. It can
  also be a character, in case it should be one of the following:
  `c("AUC", "conditionalAUC", "conditionalBoyce", "conditionalSomersD")`.

- metrics_evaluate:

  `[character]`  
  Vector of metric names to compute for model evaluation. See
  description of the argument `metric`.

- method:

  `[character="Lasso"]`  
  The penalized regression method used for fitting each model. Default
  is `method = "Lasso"`, but it could be `method = "Ridge"`,
  `"AdaptiveLasso"`, `"ElasticNet"`, or one of the ecology-constrained
  penalized regression methods (see below).

- alpha:

  `[numeric(1)=NULL]`  
  Elastic net mixing parameter. If `NULL`, glmnet chooses a default
  behavior (`alpha = 1` for Lasso, `alpha = 0` for Ridge, and
  `alpha = 0.5` for ElasticNet).

- penalty.factor:

  `[numeric,vector]`  
  Penalty factors for each coefficient in the glmnet penalty term.

- gamma:

  `[numeric(1)=1]{(0.5, 1, 2)}`  
  Gamma is the exponent for defining the vector of penalty weights when
  `method = "AdaptiveLasso`. This means that the penalties are defined
  as `penalty.factor = 1/(coef_ridge^gamma)`, where `coef_ridge` are the
  coefficients of a Ridge regression. Default is `gamma = 1`, but values
  of 0.5 or 2 could also be tried, as suggested by the authors (Zou et
  al 2006).

- standardize:

  `[character(1)="internal"]{"internal","external",FALSE}`  
  Predictor standardization mode. `"internal"` delegates standardization
  to glmnet (also standardizes dummy variables, but returns coefficients
  on original scale). `"external"` standardizes predictors before
  calling glmnet. `FALSE` skips standardization.

- predictor_table:

  `[data.frame or NULL]`  
  Default is `NULL`. Else, this is the predictor table defined by
  running `add_zoi_formula(predictor_table = TRUE)`, which is a table
  with info of all ZOI radii, shape values, and formula terms, together
  with info from other non-ZOI predictors. This table is required for
  all the ecology-constrained penalized regression methods.

- function_lasso_decay:

  `[function=log]`  
  Function used to compute decay weights for adaptive lasso penalties.
  Default is `log`.

- value_lasso_decay:

  `[numeric(1)=1]`  
  Scaling factor applied to the adaptive lasso decay function.

- function_hypothesis:

  `[function]`  
  Hypothesis function used to compute additional penalty weights.

- expected_sign_hypothesis:

  `[numeric(1)=-1]`  
  Expected coefficient sign used for hypothesis-driven penalization.

- factor_grouped_lasso:

  `[numeric(1)=1]`  
  Scaling factor applied to grouped lasso penalties.

- replace_missing_NA:

  `[logical(1)=TRUE]`  
  If `TRUE` (default), any variables missing from the data (i.e. with
  variance zero) are removed from the formula for the model fitting
  procedure, and a `NA` is set as its coefficient in the output. If
  `FALSE`, the function raises an error if there are variables with
  variance zero in the formula.

- na.action:

  `[character(1)="na.pass"]`  
  Action to take for missing values during model fitting.

- out_dir_file:

  `[character(1)=NULL]`  
  String with the prefix of the file name (and the folder) where the
  result of each model will be saved. E.g. if
  `out_dir_file = "output/test_"`, the models will be saved as RDS files
  names "test_i1.rds", "test_i2.rds", etc, within the folder "output".

- verbose:

  `[logical(1)=FALSE]`  
  Whether to print progress messages during the fitting process.

- ...:

  `[any]`  
  Options for
  [`net_logit()`](https://ninanor.github.io/oneimpact/reference/net_logit.md)
  and
  [`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html).

- phi_group:

  Additional penalty constant for the group-based penalties. A value in
  the interval 0, Inf where 0 is no additional penalty and higher values
  correspond to higher penalties.

## Value

A named list with the results for the selected metric, including:

- `coef`: coefficient matrix at optimal lambda.

- `coef_std`: standardized coefficient matrix (or `NULL`).

- `lambda`: optimal lambda value.

- `alpha`: alpha value used.

- `train_score`, `test_score`, `validation_score`: performance scores.

- `validation_score_avg`: mean validation score across blocks.

- `coefs_all`, `lambdas`: coefficients and lambdas for all evaluated
  metrics.

- `metrics_evaluated`: full detail of all evaluated metrics.

- `glmnet_fit`: the raw fitted glmnet object.

- `parms`: recorded input parameters.

A named list with the results for the selected metric, including:

- `coef`: coefficient matrix at optimal lambda.

- `coef_std`: standardized coefficient matrix (or `NULL`).

- `lambda`: optimal lambda value.

- `alpha`: alpha value used.

- `train_score`, `test_score`, `validation_score`: performance scores.

- `validation_score_avg`: mean validation score across blocks.

- `coefs_all`, `lambdas`: coefficients and lambdas for all evaluated
  metrics.

- `metrics_evaluated`: full detail of all evaluated metrics.

- `glmnet_fit`: the raw fitted glmnet object.

- `parms`: recorded input parameters.

## Details

By default, `fit_net_clogit()` does not standardize predictor variables.
If you want numeric variables to be standardized, you can either use
`[oneimpact::bag_fit_net_clogit()]` with parameter `standardize = TRUE`
or provide an already standardized data set as input.

By default, `fit_net_logit()` does not standardize predictor variables.
If you want numeric variables to be standardized, you can either use
`[oneimpact::bag_fit_net_logit()]` with parameter `standardize = TRUE`
or provide an already standardized data set as input.

## References

Zou, H., 2006. The Adaptive Lasso and Its Oracle Properties. Journal of
the American Statistical Association 101, 1418–1429.
https://doi.org/10.1198/016214506000000735

Zou, H., 2006. The Adaptive Lasso and Its Oracle Properties. Journal of
the American Statistical Association 101, 1418–1429.
https://doi.org/10.1198/016214506000000735

## See also

[`bag_fit_net_logit()`](https://ninanor.github.io/oneimpact/reference/bag_fit_net_functions.md),
[`net_logit()`](https://ninanor.github.io/oneimpact/reference/net_logit.md),
[`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html)
