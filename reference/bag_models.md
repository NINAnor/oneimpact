# Summary of a bag of models

Creates a bag of models object from a fitted ensemble of models. This
function synthesizes results from multiple model resamples, assigns
weights based on validation scores, and computes weighted coefficients
and predictions for ensemble modeling.

## Usage

``` r
bag_models(
  fitted,
  data,
  score2weight = score2weight_min_invmean,
  metric = fitted$metric,
  weights_col = c("validation_score", "habitat_validation_score")[1],
  score_threshold = 0.7,
  weights_function = NULL
)
```

## Arguments

- fitted:

  `[list]`  
  Object containing fitted models, typically output from
  `[oneimpact::fit_net_logit()]` or `[oneimpact::fit_net_clogit()]` with
  multiple resamples. Must contain elements: `$models`, `$metric`,
  `$formula`, `$method`, `$standardize`, and `$samples`.

- data:

  `[data.frame,tibble]`  
  Complete data set used to fit the models. Used here only to compute
  summary statistics and covariate scaling information.

- score2weight:

  `[function]`  
  Function to convert validation scores to model weights. Should accept
  arguments: `x` (model result), `weights_col` (column name), and
  `score_threshold`. Default is `score2weight_min_invmean`. See
  `[oneimpact::score2weight_min_invmean()]` and related functions.

- metric:

  `[character]`  
  Name of the evaluation metric to extract from models. Default uses the
  same metric from the fitted models. Common values: `"AUC"`,
  `"Cindex"`, `"conditionalAUC"`, `"conditionalSomersD"`.

- weights_col:

  `[character(1)="validation_score"]`  
  Column name containing validation scores to be converted to weights.
  One of `"validation_score"` or `"habitat_validation_score"`. The
  latter is used when habitat-specific validation scores are available,
  such as in step-selection functions in which one wants to disentangle
  the effects of movement and habitat.

- score_threshold:

  `[numeric(1)=0.7]`  
  Minimum validation score threshold. Models with validation score below
  this threshold are assigned zero weight and excluded from the bag.

- weights_function:

  `[function or NULL]`  
  Function to normalize and transform weights into the range 0,1.
  Default is `w_strech_maxmin_squared`. See
  `[oneimpact::w_strech_maxmin_squared()]` and related functions.

## Value

A `bag` object (also a `list`) containing:

- n:

  Total number of models in the bag

- n_no_errors:

  Number of successfully fitted models

- n_errors:

  Number of models that failed to fit

- n_above_threshold:

  Number of models with weights above zero

- formula:

  Model formula used for fitting

- formula_no_strata:

  Formula without strata terms

- method:

  Modeling method (e.g., "Lasso", "Ridge", "AdaptiveLasso", etc.)

- metric:

  Evaluation metric used (e.g. "AUC", "Cindex", etc.)

- weights:

  Vector of model weights (normalized)

- coef:

  Matrix of raw coefficients for each model

- coef_std:

  Matrix of standardized coefficients for each model

- wcoef:

  Weighted coefficients

- wcoef_std:

  Weighted standardized coefficients

- validation_score:

  Matrix of validation scores

- weighted_validation_score:

  Weighted validation scores

- covariate_mean_sd:

  Data frame with mean and SD of covariates

- data_summary:

  Summary statistics of variables in the original data

- errors:

  Logical vector indicating which models failed

## See also

[`fit_net_logit()`](https://ninanor.github.io/oneimpact/reference/fit_net_functions.md),
[`fit_net_clogit()`](https://ninanor.github.io/oneimpact/reference/fit_net_functions.md),
[`predict.bag()`](https://ninanor.github.io/oneimpact/reference/predict.md),
[`score2weight_min_invmean()`](https://ninanor.github.io/oneimpact/reference/bag_model_weights.md),
[`w_strech_maxmin_squared()`](https://ninanor.github.io/oneimpact/reference/bag_model_weights.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming you have fitted models from fit_net_logit or fit_net_clogit
fitted_models <- bag_fit_net_logit(f = case_ ~ x + y, data = mydata, samples = samples_list)
bag <- bag_models(fitted = fitted_models, data = mydata, metric = "AUC")
} # }
```
