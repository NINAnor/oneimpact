# Predict bag of models in space

Predict bag of models in space

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
  A bag of models, resulting from a call to
  [`bag_models()`](https://ninanor.github.io/oneimpact/reference/bag_models.md).

- data:

  `[data.frame,SpatRaster]`  
  The spatial grid or stack of rasters with the layers in space, to be
  used for the spatial prediction. All variables in the bag\$formula
  must be present in `data`.

- input_type:

  `[string(1)="df"]{"df","rast"}`  
  Type of input object. Either a `data.frame` with the predictor
  variables as columns (if `input_type = "df"`, default), or a
  `SpatRaster` with the predictor variables as layers (if
  `input_type = "rast"`). So far, only `input_type = "df"` is
  implemented.

- output_type:

  `[string(1)="rast"]{"df","rast"}`  
  Type of output object. Typically, the same type of object as
  `input_type`, but rasters can also be saved if the input is a
  `data.frame` when `output_type = "rast"`, and data.frames can also be
  saved if the input is a `SpatRaster` when `output_type = "df"`.

- prediction_type:

  `[string(1)="exp"]{"exp", "exponential", "linear"}` Type of
  transformation for the prediction. One of `"exp"` or `"expornential"`
  for exponential response or `"linear"` for linear response.

- gid:

  `[string(1)="gid"]`  
  String with the name of the "gid" or point ID column in `data`. Only
  relevant if `input_type = "df"`.

- coords:

  `[vector,string(2)=c("x", "y")]`  
  Vector with two elements with the names of the coordinates
  representing (x,y) coordinates for the pixels. Only relevant if
  `input_type = "df"`.

- crs:

  `[string(1)=NULL]`  
  Code for the coordinate reference system of the output raster. Only
  relevant if `input_type = "df"`. For more details, check
  [`terra::crs()`](https://rspatial.github.io/terra/reference/crs.html).
