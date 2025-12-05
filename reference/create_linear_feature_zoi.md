# Creates a line feature and gets the value of their zone of influence on locations over the feature

The function `create_linear_feature_zoi()` computes the cumulative zone
of influence (ZOI) of one single linear feature so it is used to
correctly create predictions and response plots for linear
infrastructure, considering the potential responses at multiple radii.

## Usage

``` r
create_linear_feature_zoi(
  radii = c(100, 250, 500, 1000, 2500, 5000, 10000),
  type = c("circle", "Gauss", "rectangle", "exp_decay", "bartlett", "threshold",
    "mfilter")[1],
  radius_max = max(radii),
  zoi_limit = 0.05,
  type_feature_recompute = FALSE,
  res = 100,
  value = 1
)
```

## Arguments

- radii:

  `[numeric,vector=c(100, 250, 500, 1000, 2500, 5000, 10000)]`  
  Vector of radii for which the zone of influence should be computed.

- type:

  `[character="circle"]{"circle", "Gauss", "rectangle", "exp_decay", "bartlett", "threshold"}`  
  Shape of the zone of influence (ZOI), Default is `circle"`. It can
  assume any of the possible values for the argument `type` in the
  function
  [`dist_decay()`](https://ninanor.github.io/oneimpact/reference/zoi_functions.md).

- radius_max:

  `[numeric=max(radii)]`  
  Maximum radius, used to set the size of the landscape/raster for ZOI
  computations.

- res:

  `[numeric(1)=100]`  
  Resolution for the raster created. This might impact what are the
  values observed in the ZOI.

- line_value:

  `[numeric(1)=1]`  
  Value set to the raster line created. Default is 1. It could be
  changed to different values if we want to represent e.g. the value in
  the linear feature as the roads traffic or another value for
  spatio-temporally dynamic variables.

## Details

This function could be extended to the nearest ZOI, if needed.

## See also

[`predict()`](https://ninanor.github.io/oneimpact/reference/predict.md)

## Examples

``` r
# create feature
create_linear_feature_zoi(radii = c(100, 250, 500, 1000, 2500, 5000, 10000),
                          type = "exp_decay",
                          res = 100)
#>       100       250       500      1000      2500      5000     10000 
#>  1.105250  1.861974  3.426254  6.690854 16.580199 33.088243 63.428301 
create_linear_feature_zoi(type = "exp_decay",
                          type_feature_recompute = FALSE)
#>       100       250       500      1000      2500      5000     10000 
#>  1.105250  1.861974  3.426254  6.690854 16.580199 33.088243 63.428301 
```
