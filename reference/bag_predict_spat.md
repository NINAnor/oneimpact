# Predict bag of models in space

Creates spatial predictions from a bag of models for a grid provided as
a data frame or raster. It can return weighted mean predictions,
weighted median/uncertainty summaries, and/or individual model
predictions.

## Usage

``` r
bag_predict_spat(
  bag,
  data,
  input_type = c("df", "rast")[1],
  output_type = c("df", "rast")[2],
  prediction_type = c("exp", "exponential", "linear")[1],
  standardize = FALSE,
  what = c("mean", "median", "ind"),
  gid = "gid",
  coords = c("x33", "y33"),
  crs = NULL,
  gridalign = FALSE,
  output_rescale = FALSE,
  prediction_max_quantile = 0.999,
  uncertainty_quantiles = c(0.25, 0.75),
  verbose = FALSE
)

bag_predict_spat_vars(
  bag,
  data,
  predictor_table_zoi,
  input_type = c("df", "rast")[1],
  output_type = c("df", "rast")[2],
  prediction_type = c("exp", "exponential", "linear")[1],
  standardize = FALSE,
  what = c("mean", "median", "ind"),
  gid = "gid",
  coords = c("x33", "y33"),
  crs = NULL,
  gridalign = FALSE,
  prediction_max_quantile = 0.999,
  uncertainty_quantiles = c(0.25, 0.75),
  verbose = FALSE
)
```

## Arguments

- bag:

  `[bag,list]`  
  A bag of models from
  [`bag_models()`](https://ninanor.github.io/oneimpact/reference/bag_models.md).

- data:

  `[data.frame,SpatRaster]`  
  Spatial grid of predictors used for prediction. All predictors in the
  model formula must be present.

- input_type:

  `[character(1)="df"]{"df","rast"}`  
  Input type: data frame with predictor values (`"df"`) or `SpatRaster`
  raster stack/layers (`"rast"`).

- output_type:

  `[character(1)="rast"]{"df","rast"}`  
  Output type: data frame or raster outputs.

- prediction_type:

  `[character(1)="exp"]{"exp","exponential","linear"}`  
  Prediction scale.

- standardize:

  `[logical(1)=FALSE]`  
  Whether to standardize predictors before prediction using bag summary
  statistics.

- what:

  `[character]`  
  Prediction outputs to compute. Any combination of "mean", "median",
  and "ind" for the weighted mean, weighted median, and individual model
  predictions, respectively.

- gid:

  `[character(1)="gid"]`  
  Name of unique pixel/grid identifier column (used for `df` input and
  output).

- coords:

  `[character(2)=c("x","y")]`  
  Names of x/y coordinate columns. Used if `input_type = "df"`.

- crs:

  `[character(1)=NULL]`  
  Coordinate reference system for raster output. Only relevant if
  `input_type = "df"`. For more details, check
  [`terra::crs()`](https://rspatial.github.io/terra/reference/crs.html).

- gridalign:

  `[logical(1)=FALSE]`  
  Whether to align output grid coordinates before rasterization.

- output_rescale:

  `[logical(1)=FALSE]`  
  Whether to rescale raster outputs to 0,1 using prediction_max_quantile
  clipping.

- prediction_max_quantile:

  `[numeric(1)=0.999]`  
  Quantile used as upper truncation threshold when
  `output_rescale = TRUE`.

- uncertainty_quantiles:

  `[numeric(2)=c(0.25, 0.75)]`  
  Quantiles used for weighted uncertainty summaries when what includes
  "median".

- verbose:

  `[logical(1)=FALSE]`  
  Print progress messages.

- predictor_table_zoi:

  `[data.frame]`  
  Predictor table with ZOI metadata (for example from
  [`add_zoi_formula()`](https://ninanor.github.io/oneimpact/reference/add_zoi_formula.md)
  predictor_table output), used to aggregate predictions by ZOI variable
  groups in `bag_predict_spat_vars()`.

## Value

For `bag_predict_spat()`, a list with prediction outputs:

- grid: prediction table in tabular form (always returned).

- weights: model weights used for prediction (non-zero weight models).

- r_weighted_avg_pred: raster of weighted mean prediction (if requested
  and `output_type = "rast"`).

- r_ind_summ_pred: raster stack with weighted median and uncertainty
  summaries, including the interquartile range and the quantile
  coefficient of variation (if requested and `output_type = "rast"`).

- r_ind_pred: raster stack of individual model predictions (if requested
  and `output_type = "rast"`).

For `bag_predict_spat_vars()`, a list with partial prediction outputs:

- vars: per-variable prediction summaries.

- grid: grid-level prediction table.

- weights: model weights used for prediction.

- r_weighted_avg_pred: raster output for weighted means (if requested).

- r_ind_summ_pred: raster output for weighted uncertainty summaries (if
  requested).

- r_ind_pred: raster output for individual-model predictions (if
  requested).
