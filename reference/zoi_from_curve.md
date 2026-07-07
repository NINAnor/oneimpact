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

## Examples
