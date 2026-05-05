# Get estimates of zone of influence (ZOI) from response curves

This generic function computes ZOI metrics (maximum effect size, ZOI
radius and impact) for ZOI predictor variables based on response curves.
The ZOI radius is estimated as the distance/radius in which the relative
selection strength corresponds to on a given percentage of the maximum
effect size (e.g. 95% ZOI radius). The impact accounts for both the
effect size and the ZOI radius and corresponds to the area under (or
over, if negative) the ZOI response curve. The function supports two
types of input: a data.frame of predictions or a bag of models.

## Usage

``` r
zoi_from_curve(x, ...)

# S3 method for class 'data.frame'
zoi_from_curve(
  x,
  weights,
  percentage = 0.95,
  curve = c("median", "mean"),
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
  curve = c("median", "mean"),
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

  Either a `data.frame` containing response curve predictions for a
  single variable, or a `bag` object containing an ensemble of models.

- ...:

  Additional arguments passed to the appropriate method.

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
  Numeric vector of quantiles used for prediction summaries.

- ci:

  `[logical(1)=TRUE]`  
  Logical. Whether to compute ZOI estimates for the upper and lower
  limits of the confidence interval. Default is TRUE.

- type:

  `[character(1)="linear"]{"linear", "exp"}`  
  Character. Defines whether the calculation of ZOI should be based on
  the prediction of at linear or response (exponential) scale:
  `"linear"` or `"exp"`, respectively.

- mean_col_name:

  `[character="mean"]`  
  Name of the column containing the mean response curve.

- median_col_name:

  `[character="quantile:0.5"]`  
  Name of the column containing the median response curve.

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
  metrics. If `TRUE`, the output is necessarily a `list` with
  predictions and the ZOI parameters.

- return_format:

  `[character="df"]{"list", "df"}`  
  Format of the returned ZOI metrics. Either a list of data.frames (if
  `return_format = "list"`), one for each variable, or a single
  `data.frame` (default, if `return_format = "df"`).

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

- ci_col_name:

  `[character=c("quantile:0.255", "quantile:0.975")]`  
  Character vector of length 2. Names of columns for lower and upper
  confidence intervals.

## Value

A `data.frame` or a `list` containing ZOI metrics:

- `max_effect_size`: Maximum effect size on the relative selection
  strength (y axis).

- `zoi_radius`: Distance at which the effect drops below a threshold,
  defined by the parameter `percentage`.

- `effect_zoi_radius`: Relative selection strength (y axis) valye in
  which the ZOI is reached.

- `impact`: Area under the curve up to the ZOI radius, combining the
  varying effect size with distance. Each ZOI measure presents mean,
  median, CI lower, and CI upper.

If `x` is a bag object, the function returns wither a `list` or
`data.frame` of ZOI measures for each ZOI variable in the bag. If
`return_predictions = TRUE`, also returns the prediction curves.

## See also

[`predict()`](https://ninanor.github.io/oneimpact/reference/predict.md),
[`plot_response()`](https://ninanor.github.io/oneimpact/reference/plot_response.md),
[`weirdness()`](https://ninanor.github.io/oneimpact/reference/weirdness.md)
\##example examples/zoi_from_curve_example.R
