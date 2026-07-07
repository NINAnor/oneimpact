# Computes ecological implausibility for a fitted model or its estimated coefficients

This function evaluates ecological plausibility in model coefficients
and response curves. Ecological plausibility refers to whether estimated
relationships between predictors and responses are consistent with prior
ecological theory, expected species–environment relationships, and
smooth asymptotic behavior. Implausible responses include abrupt sign
changes, oscillations between selection and avoidance, and coefficient
signs opposite to prior expectations.

## Usage

``` r
implausibility(x, ...)

# S3 method for class 'numeric'
implausibility(
  x,
  which_coef_sign = c("count", "sum", "raw", "index")[1],
  expected_sign = -1,
  zero_coefficient_limit = 1e-08
)

# S3 method for class 'data.frame'
implausibility(
  x,
  expected_sign = -1,
  response = c("mean", "mid")[1],
  measure = c("n_crosses", "where_crosses", "response_area_opposite",
    "response_area_ratio", "n_inflection", "difference_inflection",
    "response_area_inflection")[1]
)

# S3 method for class 'bag'
implausibility(
  x,
  data,
  measure = c("coef_sign", "n_crosses", "where_crosses", "response_area_opposite",
    "n_inflection", "difference_inflection"),
  wmean = TRUE,
  which_coef_sign = c("count", "sum")[1],
  expected_sign = -1,
  zero_coefficient_limit = 1e-08,
  which_n_cross = c("mean", "sum")[1],
  response = c("mean", "mid")[1],
  baseline = "zero",
  type_feature = "point",
  type_feature_recompute = TRUE,
  resolution = 200,
  radii = c(100, 250, 500, 1000, 2500, 5000, 10000),
  zoi_shape = c("circle", "Gauss", "rectangle", "exp_decay", "bartlett", "threshold",
    "mfilter")[1],
  radius_max = NULL,
  ...
)
```

## Arguments

- x:

  `[bag]`  
  A bag of models, resulting from a call to
  [`bag_models()`](https://ninanor.github.io/oneimpact/reference/bag_models.md).

- ...:

  `[any]`  
  Additional arguments passed to methods.

- which_coef_sign:

  `[character(1)="count"]{"count","sum","raw","index"}`  
  Which measure to use for the coefficients when
  `measure = "coef_sign"`. If `"count"` (default), only the sign matters
  and the number of coefficients with unexpected sign is returned. If
  `"sum"`, the sum of the (standardized) coefficients with unexpected
  sign is returned, accounting for magnitude.

- expected_sign:

  `[numeric(1)=-1]`  
  Expected sign of the coefficient. Either -1, +1, or 0 (no effect).

- zero_coefficient_limit:

  `[numeric(1)=1e8]`  
  Value above which an estimated coefficient is considered non-zero.
  Default is `1e-8`.

- measure:

  `[string(1)]{"coef_sign", "n_crosses", "response_area"}`  
  Measure used to quantify ecological implausibility in the model or
  coefficients. It can be one or multiple of these options:

  - `"coef_sign"`: counts coefficients whose sign is opposite to the
    ecologically expected sign.

  - `"n_crosses"`: counts sign crossings for the response curve.

  - `"response_area"`: computes area under the response curve in the
    unexpected direction.

- data:

  `[data.frame]`  
  The original, complete data used for model fitting.

## Value

The output depends on the input type and measure used. For a numeric
vector of coefficients, it returns a single value indicating the degree
of implausibility. For a data frame representing the response curves, it
returns a list with measures of ecological implausibility. For a bag of
models, it returns a list with measures of ecological implausibility for
each of the ZOI variables in the bag

## Examples

``` r
#-------
# implausibility for vector of coefficients

# implausibility for coefficients for one type of ZOI variable

# set coefficients
coefs <- c(-1, -0.5, -0.1, 0.8, 0.3, -0.1)
expected_sign <- -1
implausibility(coefs, expected_sign = expected_sign)
#> [1] 2
implausibility(coefs, expected_sign = expected_sign, which_coef = "sum")
#> [1] 1.1
implausibility(coefs, expected_sign = expected_sign, which_coef = "raw")
#> [1] 0.8 0.3
implausibility(coefs, expected_sign = expected_sign, which_coef = "index")
#> [1] 4 5

#-------
# implausibility for data.frame with (x,y) for line

# checking for lines crossing zero
x <- seq(0, 10, 0.01)
y <- -8 + 10 * x - 1.5 * x**2
df <- data.frame(x = x, y = y)
plot(x, y, ylab = "Response", xlab = "Distance from source")
abline(h = 0, col = "red")


# n crosses
implausibility(df, response = "y", measure = "n_crosses")
#> [1] 2
# where does the curve crosses zero
implausibility(df, response = "y", measure = "where_crosses")
#> [1] 0.92 5.73
# area on the opposite side of the expected sign
implausibility(df, response = "y", measure = "response_area_opposite")
#> [1] 27.7758
# ratio between area above and area on the expected sign
implausibility(df, response = "y", measure = "response_area_ratio")
#> [1] 0.2571948

# checking for inflection points
x <- seq(0, 14, 0.01)
y <- -560 + 314 * x - 56 * x**2 + 3*x**3
df <- data.frame(x = x, y = y)
plot(x, y); abline(h = 0, col = "red")

# inflection points
which(inflection(y))
#> [1] 428 819
abline(v = x[inflection(y)], lty = 2)


# n crosses
implausibility(df, response = "y", measure = "n_crosses")
#> [1] 1
# n inflection points
implausibility(df, response = "y", measure = "n_inflection")
#> [1] 2
# difference between inflection points
implausibility(df, response = "y", measure = "difference_inflection")
#> [1] 89.84515

#-------
# implausibility for bag

#---
# fit a bag to be tested

# load packages
library(glmnet)
#> Loading required package: Matrix
#> Loaded glmnet 5.0
#> 
#> Attaching package: ‘glmnet’
#> The following objects are masked from ‘package:oneimpact’:
#> 
#>     Cindex, coxnet.deviance

# load data
data("reindeer_rsf")
# rename it just for convenience
dat <- reindeer_rsf

# formula initial structure
f <- use ~ private_cabins_cumulative_XXX + public_cabins_high_cumulative_XXX +
  trails_cumulative_XXX +
  NORUTreclass +
  # poly(norway_pca_klima_axis1, 2, raw = TRUE) +
  # poly(norway_pca_klima_axis2, 2, raw = TRUE) +
  norway_pca_klima_axis1 + norway_pca_klima_axis1_sq +
  norway_pca_klima_axis2 + norway_pca_klima_axis2_sq +
  norway_pca_klima_axis3 + norway_pca_klima_axis4

# add ZOI terms to the formula
zois <- c(100, 250, 500, 1000, 2500, 5000, 10000, 20000)
ff <- add_zoi_formula(f, zoi_radius = zois, pattern = "XXX",
                      cumulative = "",
                      type = c("exp_decay"),#, "nearest_exp_decay"),
                      separator = "", predictor_table = TRUE)
f <- ff$formula
pred_table <- ff$predictor_table

# sampling - random sampling
set.seed(1234)
samples <- create_resamples(y = dat$use,
                            p = c(0.2, 0.2, 0.2),
                            times = 20,
                            colH0 = NULL)
#> [1] "Starting random sampling..."

# fit multiple models
fittedl <- bag_fit_net_logit(f,
                             data = dat,
                             samples = samples,
                             standardize = "internal", # glmnet does the standardization of covariates
                             metric = "AUC",
                             method = "AdaptiveLasso",
                             predictor_table = pred_table,
                             parallel = "mclapply",
                             mc.cores = 8) #2)

# bag models in a single object
bag_object <- bag_models(fittedl, dat, score_threshold = 0.7)

# bag_object$coef %*% bag_object$weights
# sapply(fittedl, function(x) x$train_score)

#---
# plot to check

# ZOI public cabins cumulative
dfvar = data.frame(trails_cumulative = 1e3*seq(0.2, 20, length.out = 100))

# look into curve
plot_response(bag_object,
              dfvar = dfvar,
              data = dat,
              type = "linear", zoi = TRUE,
              type_feature =  "line",
              type_feature_recompute = TRUE,
              resolution = 300,
              ci = FALSE, indiv_pred = TRUE)

# with no line, just as an example
plot_response(bag_object,
              dfvar = dfvar,
              data = dat,
              type = "linear", zoi = TRUE,
              ci = FALSE, indiv_pred = TRUE)


plot_response(bag_object,
              dfvar = dfvar,
              data = dat,
              type = "linear", zoi = TRUE,
              ci = FALSE, indiv_pred = TRUE,
              ggplot = FALSE) |>
  ggplot(aes(trails_cumulative, Resample01)) +
  geom_line()
#> Error in ggplot(plot_response(bag_object, dfvar = dfvar, data = dat, type = "linear",     zoi = TRUE, ci = FALSE, indiv_pred = TRUE, ggplot = FALSE),     aes(trails_cumulative, Resample01)): could not find function "ggplot"

# we try the function with the curve above, but then test how we could work with the more correct one below
# implausibility measures
implausibility(bag_object,
          data = dat,
          type_feature = c("point", "line", "line"))
#> Error in plot_response(x, dfvar = dfvar, data = data, type = "linear",     zoi = TRUE, type_feature_recompute = type_feature_recompute,     resolution = resolution, type_feature = type_feat, baseline = baseline,     ci = TRUE, indiv_pred = FALSE, ggplot = FALSE, ...): unused argument (ggplot = FALSE)

# for each individual model
implausibility(bag_object,
          data = dat,
          wmean = FALSE)
#> Error in plot_response(x, dfvar = dfvar, data = data, type = "linear",     zoi = TRUE, type_feature_recompute = type_feature_recompute,     resolution = resolution, type_feature = type_feat, baseline = baseline,     wq_probs = NULL, ci = FALSE, indiv_pred = TRUE, ggplot = FALSE,     ...): unused argument (ggplot = FALSE)
```
