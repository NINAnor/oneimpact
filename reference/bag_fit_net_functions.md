# Fit a bag of conditional logistic regression/SSF/iSSF models with penalized regression in a train-validate-test setup

Fits a bag on conditional logistic regressions by calling
[`fit_net_clogit()`](https://ninanor.github.io/oneimpact/reference/fit_net_functions.md)
multiple times to fit different resamples of the input data.

Fits a bag of logistic regressions by calling
[`fit_net_logit()`](https://ninanor.github.io/oneimpact/reference/fit_net_functions.md)
multiple times to fit different resamples of the input data.

## Usage

``` r
bag_fit_net_clogit(
  f,
  data,
  samples,
  subset_samples = 1:length(samples$train),
  kernel_vars = c("step_length", "ta"),
  metric = c("coxnet.deviance", "Cindex", "conditionalAUC", "conditionalSomersD")[1],
  metrics_evaluate = c("coxnet.deviance", "Cindex", "conditionalAUC"),
  method = c("Lasso", "Ridge", "AdaptiveLasso", "DistanceDecayLasso", "DDLasso",
    "OneZOI-AdaptiveLasso", "OZ-AdaptiveLasso", "Grouped-AdaptiveLasso",
    "G-AdaptiveLasso", "ElasticNet")[1],
  standardize = c("internal", "external", FALSE)[1],
  alpha = NULL,
  penalty.factor = NULL,
  predictor_table = NULL,
  na.action = "na.pass",
  out_dir_file = NULL,
  parallel = c(FALSE, "foreach", "mclapply")[1],
  mc.cores = 2L,
  verbose = FALSE,
  ...
)

bag_fit_net_logit(
  f,
  data,
  samples,
  subset_samples = 1:length(samples$train),
  metric = c(AUC, conditionalBoyce, conditionalSomersD, conditionalAUC)[[1]],
  metrics_evaluate = c("AUC"),
  method = c("Lasso", "Ridge", "AdaptiveLasso", "DistanceDecayLasso", "DDLasso",
    "TruncatedLasso", "OneZOI-AdaptiveLasso", "OZ-AdaptiveLasso",
    "Grouped-AdaptiveLasso", "G-AdaptiveLasso", "ElasticNet")[1],
  standardize = c("internal", "external", FALSE)[1],
  alpha = NULL,
  penalty.factor = NULL,
  predictor_table = NULL,
  na.action = "na.pass",
  out_dir_file = NULL,
  parallel = c(FALSE, "foreach", "mclapply")[1],
  mc.cores = 2L,
  verbose = FALSE,
  ...
)
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
  validate. Typically computed by
  [`create_resamples()`](https://ninanor.github.io/oneimpact/reference/create_resamples.md).

- subset_samples:

  `[vector]`  
  Vector of sample indices to run (e.g. `c(1,2,4,5)` or `3:10`). By
  default, all samples in `samples` are run.

- kernel_vars:

  `[character]`  
  Names of movement kernel variables in the model (e.g. `"step_length"`
  and `"turning_angle"`).

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

  `[character(1)="Lasso"]`  
  Penalized regression method. Default is `method = "Lasso"`, but it
  could be `method = "Ridge"`, `"AdaptiveLasso"`, `"ElasticNet"`, or
  different versions of ecology-informed methods. See
  [`fit_net_logit()`](https://ninanor.github.io/oneimpact/reference/fit_net_functions.md)
  for valid options.

- standardize:

  `[character(1)="internal"]{"internal","external",FALSE}`  
  Predictor standardization mode. `"internal"` uses glmnet's internal
  standardization; `"external"` standardizes before calling glmnet;
  `FALSE` skips standardization, in case which the standardization
  should have been done before calling the function to fit the models.

- alpha:

  `[numeric(1)=NULL]`  
  Elastic net mixing parameter. If `NULL`, set automatically by method.

- penalty.factor:

  `[numeric]`  
  Penalty factors per coefficient. If `NULL`, set automatically by
  method.

- predictor_table:

  `[data.frame,NULL]`  
  Predictor table with ZOI metadata, required for distance-decay and
  ecology-constrained methods. Typically created through
  `[oneimpact::add_zoi_formula()]`.

- na.action:

  `[character(1)="na.pass"]`  
  Action to take for missing values during fitting.

- out_dir_file:

  `[character(1)=NULL]`  
  Prefix path for saving individual model results as RDS files.

- parallel:

  `[character(1)=FALSE]{"FALSE","foreach","mclapply"}`  
  Parallelization backend. If FALSE (default), no parallelization is
  made.

- mc.cores:

  `[integer(1)=2]`  
  Number of cores for `mclapply`. If `parallel == "foreach"`, cores must
  be registered beforehand using
  [`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html)
  and `doParallel::registerDoParallel()`.

- verbose:

  `[logical(1)=FALSE]`  
  Whether to print progress messages during fitting.

- ...:

  `[any]`  
  Additional options passed to
  [`fit_net_logit()`](https://ninanor.github.io/oneimpact/reference/fit_net_functions.md),
  [`net_logit()`](https://ninanor.github.io/oneimpact/reference/net_logit.md)
  and
  [`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html).

## Value

A list with:

- `n`: total number of models (resamples).

- `formula`: model formula.

- `method`: fitting method used.

- `metric`: metric used for lambda selection.

- `samples`: resampling structure.

- `standardize`: standardization mode used.

- `models`: named list of fitted model objects, one per resample.

- `covariate_mean_sd`: data frame of predictor means and SDs used for
  external standardization (or `NULL`).

- `numeric_covs`: logical vector identifying which covariates are
  numeric.

A list with:

- `n`: total number of models (resamples).

- `formula`: model formula.

- `method`: fitting method used.

- `metric`: metric used for lambda selection.

- `samples`: resampling structure.

- `standardize`: standardization mode used.

- `models`: named list of fitted model objects, one per resample.

- `covariate_mean_sd`: data frame of predictor means and SDs used for
  external standardization (or `NULL`).

- `numeric_covs`: logical vector identifying which covariates are
  numeric.

## See also

[`fit_net_clogit()`](https://ninanor.github.io/oneimpact/reference/fit_net_functions.md),
[`net_clogit()`](https://ninanor.github.io/oneimpact/reference/net_clogit.md),
[`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html)

[`fit_net_logit()`](https://ninanor.github.io/oneimpact/reference/fit_net_functions.md),
[`net_logit()`](https://ninanor.github.io/oneimpact/reference/net_logit.md),
[`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html)
