# Fits a conditional logistic regression/SSF/iSSF with penalized regression using glmnet in a train-validate-test setup

By default, `fit_net_clogit()` does not standardize predictor variables.
If you want numeric variables to be standardized, you can either use
`[bag_fit_net_clogit()]` with parameter `standardize = TRUE` or provide
an already standardized data set as input.

By default, `fit_net_logit()` does not standardize predictor variables.
If you want numeric variables to be standardized, you can either use
`[bag_fit_net_logit()]` with parameter `standardize = TRUE` or provide
an already standardized data set as input.

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
    "G-AdaptiveLasso", "HypothesisDriven-AdaptiveLasso", "HD-AdaptiveLasso",
    "ElasticNet")[1],
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
    "G-AdaptiveLasso", "HypothesisDriven-AdaptiveLasso", "HD-AdaptiveLasso",
    "ElasticNet")[1],
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
    "G-AdaptiveLasso", "HypothesisDriven-AdaptiveLasso", "HD-AdaptiveLasso",
    "ElasticNet")[1],
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
    "Grouped-AdaptiveLasso", "G-AdaptiveLasso", "HypothesisDriven-AdaptiveLasso",
    "HD-AdaptiveLasso", "ElasticNet")[1],
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
    "Grouped-AdaptiveLasso", "G-AdaptiveLasso", "HypothesisDriven-AdaptiveLasso",
    "HD-AdaptiveLasso", "ElasticNet")[1],
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

- method:

  `[character="Lasso"]`  
  The penalized regression method used for fitting each model. Default
  is `method = "Lasso"`, but it could be `method = "Ridge"` or different
  flavors of `"AdaptiveLasso"` (see details below).

- gamma:

  `[numeric(1)=1]{(0.5, 1, 2)}`  
  Gamma is the exponent for defining the vector of penalty weights when
  `method = "AdaptiveLasso`. This means that the penalties are defined
  as `penalty.factor = 1/(coef_ridge^gamma)`, where `coef_ridge` are the
  coefficients of a Ridge regression. Default is `gamma = 1`, but values
  of 0.5 or 2 could also be tried, as suggested by the authors (Zou et
  al 2006).

- standardize:

  `[logical(1)=TRUE]`  
  Logical flag for predictor variable standardization, prior to fitting
  the model sequence. The coefficients are always returned on the
  original scale. Default is standardize=TRUE. If variables are in the
  same units already, you might not wish to standardize them.

- replace_missing_NA:

  `[logical(1)=TRUE]`  
  If `TRUE` (default), any variables missing from the data (i.e. with
  variance zero) are removed from the formula for the model fitting
  procedure, and a `NA` is set as its coefficient in the output. If
  `FALSE`, the function raises an error if there are variables with
  variance zero in the formula.

- out_dir_file:

  `[character(1)=NULL]`  
  String with the prefix of the file name (and the folder) where the
  result of each model will be saved. E.g. if
  `out_dir_file = "output/test_"`, the models will be saved as RDS files
  names "test_i1.rds", "test_i2.rds", etc, within the folder "output".

- ...:

  Options for
  [`net_logit()`](https://ninanor.github.io/oneimpact/reference/net_logit.md)
  and
  [`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html).

- phi_group:

  Additional penalty constant for the group-based penalties. A value in
  the interval 0, Inf where 0 is no additional penalty and higher values
  correspond to higher penalties.

## References

Zou, H., 2006. The Adaptive Lasso and Its Oracle Properties. Journal of
the American Statistical Association 101, 1418–1429.
https://doi.org/10.1198/016214506000000735

Zou, H., 2006. The Adaptive Lasso and Its Oracle Properties. Journal of
the American Statistical Association 101, 1418–1429.
https://doi.org/10.1198/016214506000000735
