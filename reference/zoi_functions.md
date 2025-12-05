# Zone of Influence (ZOI) functions

Computes the decay functions that represent the Zone of Influence (ZOI).
The functions' radius (parameter `radius`) controls how far the zone of
influence of an infrastructure or disturbance reaches, and the
functions' shape (parameter `type`) represent represent how the ZOI
decays in space. Given a function shape (`type`) is chosen, the rate of
decay of the different ZOI functions is parameterized based on the ZOI
radius – e.g the slope of `linear_decay()` is defined so that the
function decreases to zero at the ZOI radius. These functions can be
used to transform arrays of (Euclidean) distance values (in one
dimension) or rasters of (Euclidean) distance (in two dimensions) into
ZOI values. The distances might represent the distance to human
infrastructure, sources of disturbance, or more broadly any type of land
use class or spatial variable.

## Usage

``` r
dist_decay(
  x,
  radius = NULL,
  type = c("exp_decay", "gaussian_decay", "linear_decay", "threshold_decay")[1],
  zoi_limit = 0.05,
  origin = 0,
  oneside = TRUE,
  ...
)

threshold_decay(x, radius, intercept = 1, origin = 0, oneside = TRUE)

# S3 method for class 'numeric'
threshold_decay(x, radius, intercept = 1, origin = 0, oneside = TRUE)

# S3 method for class 'SpatRaster'
threshold_decay(x, radius, intercept = 1, origin = 0, oneside = TRUE)

step_decay(x, radius, intercept = 1, origin = 0, oneside = TRUE)

bartlett_decay(x, radius, intercept = 1, origin = 0, oneside = TRUE)

# S3 method for class 'numeric'
bartlett_decay(x, radius, intercept = 1, origin = 0, oneside = TRUE)

# S3 method for class 'SpatRaster'
bartlett_decay(x, radius, intercept = 1, origin = 0, oneside = TRUE)

tent_decay(x, radius, intercept = 1, origin = 0, oneside = TRUE)

linear_decay(x, radius, intercept = 1, origin = 0, oneside = TRUE)

gaussian_decay(
  x,
  radius = NULL,
  zoi_limit = 0.05,
  intercept = 1,
  lambda = NULL,
  sigma = NULL,
  origin = 0,
  ...
)

half_norm_decay(
  x,
  radius = NULL,
  zoi_limit = 0.05,
  intercept = 1,
  lambda = NULL,
  sigma = NULL,
  origin = 0,
  ...
)

exp_decay(
  x,
  radius = NULL,
  zoi_limit = 0.05,
  intercept = 1,
  lambda = NULL,
  origin = 0,
  oneside = TRUE,
  half_life = NULL,
  zoi_hl_ratio = NULL
)
```

## Arguments

- x:

  `[numeric,SpatRaster,RasterLayer]`  
  Euclidean distance from an infrastructure, source of disturbance, or
  feature/class of interest. It can be a single value, an array of
  values, or a raster object. It must not necessarily be an Euclidean
  distance, but preferably it should be a distance measured in meters,
  to ease interpretation (e.g. geodesic distance).

- radius:

  `[numeric(1)]`  
  Radius of the zone of influence (ZOI), the distance at which the ZOI
  vanishes or goes below a given minimum limit value `zoi_limit`. See
  details.

- type:

  `[character(1)="Gauss"]{"Gauss", "exp_decay", "bartlett", "linear", "tent", "threshold", "step"}`  
  Type or shape of the decay distance.

  - If `"Gauss"` or `"half_norm"`, the ZOI follows a half-normal
    shape:  
    `intercept * exp(-lambda * (euclidean_distance^2))`. `intercept` and
    `lambda` are parameters to be defined – see details.

  - If `"exp_decay"`, the ZOI follows an exponential decay shape:  
    `intercept * exp(-lambda * euclidean_distance)`. `intercept` and
    `lambda` are parameters to be defined – see details.

  - If `"bartlett"`, `"linear_decay"`, or `"tent_decay"`, the ZOI
    follows a linear decay shape within the ZOI radius (parameter
    `radius`).

  - If `"threshold"` or `"step"`, a constant influence is consider
    within the zone of influence radius (parameter `radius`). All pixels
    closer than `radius` to infrastructure are considered as "under the
    influence" of the nearest feature, with a constant influence value
    defined by the `intercept` parameter, and all values/pixels beyond
    `radius` are assumed to have zero influence.

- zoi_limit:

  `[numeric(1)=0.05]`  
  For non-vanishing functions (e.g. `exp_decay`, `gaussian_decay`), this
  value is used to set the relationship between the ZOI radius and the
  decay functions: `radius` is defined as the minimum distance `x` at
  which the ZOI assumes values below `zoi_limit`. The default is 0.05.
  This parameter is used only if `radius` is not `NULL`.

- origin:

  `[numeric(1)=0]`  
  In which position (in 1 dimension) is located the infrastructure or
  source of disturbance? Default is zero. For raster objects, this
  parameter should be ignored.

- oneside:

  `[logical(1)=TRUE]`  
  If `FALSE`, negative distance values are considered symmetrically and
  their transformation is always positive. This parameter is only
  meaningful if `x` is a vector of values, not a raster object.

- intercept:

  `[numeric(1)=1]`  
  Maximum value of the ZOI function at when the distance from
  disturbance sources is zero (`x = 0`). For the `threshold_decay` and
  `step_decay` functions, `intercept` is the constant value of the Zone
  of Influence within the ZOI `radius`. For the other ZOI functions,
  `intercept` is the value of the functions at the origin (where the
  sources of disturbance are located, i.e. `x = 0`). Default is
  `intercept = 1`.

- lambda:

  `[numeric(2)=NULL]`  
  For the `gaussian_decay` and `exp_decay` functions, `lambda` is the
  decay parameter of the Gaussian or exponential decay function. Notice
  that the interpretation of `lambda` is different depending on the the
  function – see details for definitions. For the Gaussian decay
  function, the value for `lambda` is only considered if both
  `radius = NULL` and `sigma = NULL`. For the exponential decay
  function, the value for `lambda` is only considered if both
  `radius = NULL` and `half_life = NULL`.

- sigma:

  `[numeric(1)=NULL]`  
  Standard deviation of the Gaussian decay function. It is related to
  the Gaussian decay rate \\\lambda\\ as `lambda = 1/(2*sigma^2)`. Only
  considered to compute the ZOI for the `gaussian_decay` function when
  the ZOI radius parameter is null (`radius = NULL`).

- half_life:

  `[numeric(1)=NULL]`  
  Half life of the exponential decay function, in meters (or map units,
  for rasters). By definition, the half life is the distance where the
  exponential decay function reaches 0.5 of its maximum value. For the
  `exp_decay` function, if the ZOI radius parameter is null
  (`radius = NULL`), the value of the exponential half life
  (`half_life = log(2)/lambda`) can be used to parameterize the
  exponential decay function.

- zoi_hl_ratio:

  `[numeric(1)=NULL]`  
  For the `exp_decay` function, if both the ZOI radius `radius` and
  `zoi_hl_ratio` are given and `half_life` is `NULL`, this value is used
  to set the ZOI radius (and `zoi_limit` is ignored). `zoi_hl_ratio` is
  the ratio between the ZOI radius value and the half life of the
  exponential function. For instance, if `radius = 1200` and
  `zoi_hl_ratio = 6`, this means `half_life` is 200. As a consequence,
  the exponential decay ZOI function decreases to 0.5 at distance 200,
  and the ZOI radius = 1200 is defined as the distance at which the ZOI
  decreases to 0.5\*\*6 = 0.015625.

## Value

The ZOI values for a given array of distance values if `x` is numeric,
or a raster object delimiting the ZOI if `x` corresponds to the distance
from infrastructure or disturbance sources in 2 dimensions space.

## Details

A generic function `dist_decay()` can be used to compute ZOI values
according to functions with different shapes (parameter `type`) and
radii (parameter `radius`). Alternatively, there are specific functions
implemented for each ZOI shape.

For the threshold function (`threshold_decay()`) and the linear decay
function (`linear_decay()`), the ZOI radius (parameter `radius`) is the
distance `x` where the ZOI function value decreases to zero. For the
linear decay, this is done by setting the slope of the linear function
as `-intercept/radius`, where `intercept` is the intercept of the linear
function (here, the maximum value at `x = 0`).

For non-vanishing functions that approach zero asymptotically
(`exp_decay()`, `gaussian_decay()`), a certain limit value must be given
to define the ZOI radius – so that the ZOI radius is defined as the
distance `x` where the ZOI function goes below this limit value. For
these functions, different parameters are available for setting the
relationship between the ZOI function value and the ZOI radius.

Some functions have multiple possible names, for the sake of
flexibility:

- `linear_decay()`, `bartlett_decay()`, and `tent_decay()` are the same
  function;

- `threshold_decay()` and `step_decay()` are the same function;

- `gaussian_decay()` and `half_norm_decay()` are the same function.

Alternatively, `dist_decay()` can call all of them, given a ZOI shape is
specified through the parameter `type`.

Other functions might be implemented in the future.

## Definitions

Here are some formal definitions for the ZOI functions \\\phi(d_i, r)\\,
where \\d_i\\ is the distance to the feature \\i\\ of an infrastructure
or source of disturbance and \\r\\ is the radius of the zone of
influence:

- `threshold_decay()`: the threshold or step decay function
  \\\phi\_{threshold}\\ is positive and constant within the ZOI radius
  \\r\\, and null for \\x \ge r\\: \$\$ \phi\_{threshold}(d_i, r_k) = c
  \text{ if } d_i \< r, 0 \text{ otherwise} \$\$ where \\c\\ is a
  constant value (by default `c = 1`).

- `linear_decay()`: the linear (or tent/Bartlett) decay function
  \\\phi\_{linear}\\ decreases linearly from a maximum value \\c\\ (the
  intercept, by default `c = 1`) to zero when \\x \ge r\\:
  \$\$\phi\_{linear}(d_i, r) = c - c/r \text{ if } x \< r, \text{ 0
  otherwise}\$\$

- `exp_decay()`: the exponential decay function \\\phi\_{exp}\\
  decreases exponentially from a maximum value \\c\\ (by default
  `c = 1`) with a rate \\\lambda\\, which is defined by \\r\\ and a ZOI
  limit value \\\phi\_{lim}\\, a small ZOI value below which the
  influence is considered negligible: \$\$\phi\_{exp}(d_i, r,
  \phi\_{lim}) = c exp(-\lambda d_i)\$\$ with \$\$\lambda =
  ln(1/\phi\_{lim}) / r\$\$ In this context, the ZOI radius \\r\\ is the
  distance beyond which \\\phi\_{exp} \< \phi\_{lim}\\.

- `gaussian_decay()`: the Gaussian decay function \\\phi\_{Gauss}\\
  follows a Gaussian (half-normal) decay with maximum \\c\\ (by default
  `c = 1`) and a decay rate \\\lambda\\ defined by \\r\\ and a ZOI limit
  value \\\phi\_{lim}\\, a small ZOI value below which the influence is
  considered negligible: \$\$\phi\_{Gauss}(d_i, r, \phi\_{lim}) = c
  exp(-\lambda d_i^2)\$\$ with \$\$\lambda = ln(1/\phi\_{lim}) /
  (r^2)\$\$ In this context, the ZOI radius \\r\\ is the distance beyond
  which \\\phi\_{exp} \< \phi\_{lim}\\. Note that \\\lambda\\ is defined
  differently for the `gaussian_decay` and the `exp_decay` functions.

## Parameterization

Some of the shapes of the ZOI (parameter `type` in `dist_decay()`) might
be parameterized in multiple ways. Here is a brief description of each
possibility:

- For the `"Gauss"` or `"half_norm"` shapes, the ZOI follows a
  half-normal shape:  
  `intercept * exp(-lambda * (euclidean_distance^2))`. `intercept` and
  `lambda` are parameters to be defined. There are three ways of
  specifying `lambda`:

  - If the `radius = NULL` (default), `lambda` is a parameter by itself
    to be specified by the user. In all other cases (below) the value of
    this parameter is ignored, even if provided.

  - If the parameter `radius` is provided, the rate of decay is given by
    `lambda = log(1/zoi_limit) / (radius**2)`. In other words, `lambda`
    is defined so that the function decreases to `zoi_limit` when
    `x = radius`.

  - If the `radius = NULL` and `sigma` is provided, `lambda` is defined
    as `lambda = 1/(2*sigma**2)`.

- For the `"exp_decay"` shape, the ZOI follows an exponential decay
  shape:  
  `intercept * exp(-lambda * euclidean_distance)`. `intercept` and
  `lambda` are parameters to be defined. There are four ways of
  specifying `lambda`:

  - If the `radius = NULL` (default), `lambda` is a parameter by itself
    to be specified by the user. In all other cases (below) the value of
    this parameter is ignored, even if provided.

  - If the parameter `radius` is provided, the rate of decay is given by
    `lambda = log(1/zoi_limit) / radius`. In other words, `lambda` is
    defined so that the function decreases to `zoi_limit` when
    `x = radius`.

  - If the `radius = NULL` and `half_life` is given, `lambda` is defined
    based on the half life of the exponential function – the distance at
    which the function decreases to 1/2. If `zoi_hl_ratio = NULL`,
    `lambda` is defined as `lambda = log(2)/half_life`.

  - The last possibility is to specify `zoi_hl_ratio`, the ratio between
    the ZOI radius and the half life of the exponential function. For
    instance, if `zoi_hl_ratio = 4`, this means the ZOI radius is
    defined as `4*half_life`. If `zoi_hl_ratio` is provided, the
    exponential `half_life` is defined based on this parameter and
    `lambda` is defined accordingly, based on the relationship above. In
    this case, `radius` is ignored, even if specified.

- For the `"bartlett"`, `"linear_decay"`, or `"tent_decay"` shapes, the
  ZOI follows a linear decay shape (`y = a*x + b`) within the ZOI radius
  (parameter `radius`). The intercept of the linear function (`b`) is
  given by the parameter `intercept` and the slope (`a`) is given by
  `-intercept/radius`.

- For the `"threshold"` or `"step"` shapes, a constant influence is
  consider within the zone of influence radius (parameter `radius`). All
  pixels closer than `radius` to infrastructure are considered as "under
  the influence" of the nearest feature, with a constant influence value
  defined by the `intercept` parameter, and all values/pixels beyond
  `radius` are assumed to have zero influence.

## Examples

``` r
# generic dist_decay function
oneimpact::dist_decay(500, radius = 1000, type = "exp_decay")
#> [1] 0.2236068
oneimpact::dist_decay(500, radius = 1000, type = "gaussian_decay")
#> [1] 0.4728708
oneimpact::dist_decay(500, radius = 1000, type = "linear_decay")
#> [1] 0.5
oneimpact::dist_decay(500, radius = 1000, type = "step_decay")
#> [1] 1

# test the zone of influence functions
# here we use ggplot() to illustrate the functions, to make the figures more
# widely reproducible
# to ease the plots, use the function plot_zoi1d()
library(ggplot2)

# exponential decay
exp_decay(10, radius = 30)
#> [1] 0.3684031

f1 <- ggplot(data.frame(x = c(0, 30)), aes(x = x)) +
  stat_function(fun = exp_decay, args = list(radius = 20)) +
  labs(x = "Distance", y = "Zone of Influence") +
  theme_bw()
f1


# exponential decay - two sides
f1_2 <- ggplot(data.frame(x = c(-30, 30)), aes(x = x)) +
  stat_function(fun = exp_decay,
                args = list(radius = 20, oneside = FALSE)) +
  labs(x = "Distance", y = "Zone of Influence") +
  theme_bw()
f1_2


# threshold
threshold_decay(5, radius = 10)
#> [1] 1
threshold_decay(10, radius = 10)
#> [1] 0

f2 <- ggplot(data.frame(x = c(0, 30)), aes(x = x)) +
  stat_function(fun = threshold_decay,
                args = list(radius = 20), linetype = 2) +
  labs(x = "Distance", y = "Zone of Influence") +
  theme_bw()
f2


# threshold - two sides
f2_2 <- ggplot(data.frame(x = c(-30, 50)), aes(x = x)) +
  stat_function(fun = threshold_decay,
                args = list(radius = 20, oneside = FALSE), linetype = 2) +
  labs(x = "Distance", y = "Zone of Influence") +
  theme_bw()
f2_2


# linear, tent, or bartlett decay
bartlett_decay(5, radius = 10)
#> [1] 0.5
bartlett_decay(8, radius = 10)
#> [1] 0.2

f3 <- ggplot(data.frame(x = c(0, 30)), aes(x = x)) +
  stat_function(fun = bartlett_decay, args = list(radius = 20), linetype = 3) +
  labs(x = "Distance", y = "Zone of Influence") +
  theme_bw()
f3


# linear, two sides
f3_3 <- ggplot(data.frame(x = c(-30, 40)), aes(x = x)) +
  stat_function(fun = bartlett_decay,
                args = list(radius = 20, origin = 10, oneside = FALSE), linetype = 3) +
  labs(x = "Distance", y = "Zone of Influence") +
  theme_bw()
f3_3


# guassian or half normal
gaussian_decay(5, sigma = 6)
#> [1] 0.7066483

f4 <- ggplot(data.frame(x = c(0, 30)), aes(x = x)) +
  stat_function(fun = gaussian_decay,
                args = list(radius = 20, zoi_limit = 0.05), linetype = 4) +
  labs(x = "Distance", y = "Zone of Influence") +
  geom_vline(xintercept = 20, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0.05, linetype = 2, color = "darkgrey") +
  theme_bw()
f4


# half normal - two sides
gaussian_decay(5, sigma = 6)
#> [1] 0.7066483

f4_2 <- ggplot(data.frame(x = c(-30, 30)), aes(x = x)) +
  stat_function(fun = gaussian_decay,
                args = list(radius = 20, zoi_limit = 0.05), linetype = 4) +
  labs(x = "Distance", y = "Zone of Influence") +
  geom_vline(xintercept = c(-20, 20), linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0.05, linetype = 2, color = "darkgrey") +
  theme_bw()
f4_2


# plot several ZoI with the same radius
f1 +
  stat_function(fun = threshold_decay, args = list(radius = 20), linetype = 2) +
  stat_function(fun = bartlett_decay, args = list(radius = 20), linetype = 3) +
  stat_function(fun = gaussian_decay, args = list(radius = 20, zoi_limit = 0.05), linetype = 4) +
  labs(x = "Distance", y = "Zone of Influence") +
  theme_bw()


#---
# applying dist_decay functions for rasters
library(terra)

# calculate Euclidean distance
f <- system.file("raster/sample_area_cabins.tif", package = "oneimpact")
cabins <- terra::rast(f)
cabins_dist <- calc_zoi_nearest(cabins, type = "euclidean")

# transform Euclidean in distance decay
# exponential decay
plot(oneimpact::dist_decay(cabins_dist, radius = 1000, type = "exp_decay"))

# linear decay
plot(oneimpact::dist_decay(cabins_dist, radius = 1000, type = "tent_decay"))
```
