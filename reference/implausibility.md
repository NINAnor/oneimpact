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
#> $n_coefs
#> [1] 24
#> 
#> $n_resamples
#> [1] 19
#> 
#> $coef_sign_index
#> $coef_sign_index$private_cabins_cumulative_
#> $coef_sign_index$private_cabins_cumulative_[[1]]
#> [1] 5 7
#> 
#> 
#> $coef_sign_index$public_cabins_high_cumulative_
#> $coef_sign_index$public_cabins_high_cumulative_[[1]]
#> [1] 3 7
#> 
#> 
#> $coef_sign_index$trails_cumulative_
#> $coef_sign_index$trails_cumulative_[[1]]
#> [1] 2 6
#> 
#> 
#> 
#> $coef_sign_names
#> $coef_sign_names$private_cabins_cumulative_
#> $coef_sign_names$private_cabins_cumulative_[[1]]
#> [1] "private_cabins_cumulative_exp_decay2500" 
#> [2] "private_cabins_cumulative_exp_decay10000"
#> 
#> 
#> $coef_sign_names$public_cabins_high_cumulative_
#> $coef_sign_names$public_cabins_high_cumulative_[[1]]
#> [1] "public_cabins_high_cumulative_exp_decay500"  
#> [2] "public_cabins_high_cumulative_exp_decay10000"
#> 
#> 
#> $coef_sign_names$trails_cumulative_
#> $coef_sign_names$trails_cumulative_[[1]]
#> [1] "trails_cumulative_exp_decay250"  "trails_cumulative_exp_decay5000"
#> 
#> 
#> 
#> $coef_sign_radii
#> $coef_sign_radii$private_cabins_cumulative_
#> $coef_sign_radii$private_cabins_cumulative_[[1]]
#> [1]  2500 10000
#> 
#> 
#> $coef_sign_radii$public_cabins_high_cumulative_
#> $coef_sign_radii$public_cabins_high_cumulative_[[1]]
#> [1]   500 10000
#> 
#> 
#> $coef_sign_radii$trails_cumulative_
#> $coef_sign_radii$trails_cumulative_[[1]]
#> [1]  250 5000
#> 
#> 
#> 
#> $coef_sign_value
#> $coef_sign_value$private_cabins_cumulative_
#> $coef_sign_value$private_cabins_cumulative_[[1]]
#> [1] 0.077101490 0.005901955
#> 
#> 
#> $coef_sign_value$public_cabins_high_cumulative_
#> $coef_sign_value$public_cabins_high_cumulative_[[1]]
#> [1] 0.2458388 3.2986788
#> 
#> 
#> $coef_sign_value$trails_cumulative_
#> $coef_sign_value$trails_cumulative_[[1]]
#> [1] 0.007444606 0.060576073
#> 
#> 
#> 
#> $coef_sign
#>     private_cabins_cumulative_ public_cabins_high_cumulative_ 
#>                              2                              2 
#>             trails_cumulative_ 
#>                              2 
#> 
#> $coef_sign_sum
#> [1] 6
#> 
#> $cross_index
#> $cross_index$private_cabins_cumulative_
#> integer(0)
#> 
#> $cross_index$public_cabins_high_cumulative_
#> integer(0)
#> 
#> $cross_index$trails_cumulative_
#> integer(0)
#> 
#> 
#> $where_crosses
#> $where_crosses$private_cabins_cumulative_
#> list()
#> 
#> $where_crosses$public_cabins_high_cumulative_
#> list()
#> 
#> $where_crosses$trails_cumulative_
#> list()
#> 
#> 
#> $n_crosses
#>     private_cabins_cumulative_ public_cabins_high_cumulative_ 
#>                              0                              0 
#>             trails_cumulative_ 
#>                              0 
#> 
#> $n_crosses_total
#> [1] 0
#> 
#> $response_area_opposite
#>     private_cabins_cumulative_ public_cabins_high_cumulative_ 
#>                              0                              0 
#>             trails_cumulative_ 
#>                              0 
#> 
#> $response_area_opposite_total
#> [1] 0
#> 
#> $response_area_ratio
#>     private_cabins_cumulative_ public_cabins_high_cumulative_ 
#>                              0                              0 
#>             trails_cumulative_ 
#>                              0 
#> 
#> $response_area_ratio_total
#> [1] 0
#> 
#> $n_inflection
#>     private_cabins_cumulative_ public_cabins_high_cumulative_ 
#>                              0                              0 
#>             trails_cumulative_ 
#>                              2 
#> 
#> $n_inflection_total
#> [1] 2
#> 
#> $difference_inflection
#>     private_cabins_cumulative_ public_cabins_high_cumulative_ 
#>                      0.0000000                      0.0000000 
#>             trails_cumulative_ 
#>                      0.1090383 
#> 
#> $difference_inflection_total
#> [1] 0.1090383
#> 

# for each individual model
implausibility(bag_object,
          data = dat,
          wmean = FALSE)
#> $n_coefs
#> [1] 24
#> 
#> $n_resamples
#> [1] 19
#> 
#> $coef_sign_index
#> $coef_sign_index$private_cabins_cumulative_
#> $coef_sign_index$private_cabins_cumulative_[[1]]
#> integer(0)
#> 
#> $coef_sign_index$private_cabins_cumulative_[[2]]
#> [1] 7
#> 
#> $coef_sign_index$private_cabins_cumulative_[[3]]
#> integer(0)
#> 
#> $coef_sign_index$private_cabins_cumulative_[[4]]
#> [1] 5 7
#> 
#> $coef_sign_index$private_cabins_cumulative_[[5]]
#> [1] 5 7
#> 
#> $coef_sign_index$private_cabins_cumulative_[[6]]
#> [1] 7
#> 
#> $coef_sign_index$private_cabins_cumulative_[[7]]
#> integer(0)
#> 
#> $coef_sign_index$private_cabins_cumulative_[[8]]
#> [1] 5 7
#> 
#> $coef_sign_index$private_cabins_cumulative_[[9]]
#> [1] 5 7
#> 
#> $coef_sign_index$private_cabins_cumulative_[[10]]
#> [1] 7
#> 
#> $coef_sign_index$private_cabins_cumulative_[[11]]
#> [1] 5 7
#> 
#> $coef_sign_index$private_cabins_cumulative_[[12]]
#> integer(0)
#> 
#> $coef_sign_index$private_cabins_cumulative_[[13]]
#> [1] 5 7
#> 
#> $coef_sign_index$private_cabins_cumulative_[[14]]
#> [1] 7
#> 
#> $coef_sign_index$private_cabins_cumulative_[[15]]
#> [1] 7
#> 
#> $coef_sign_index$private_cabins_cumulative_[[16]]
#> [1] 5 7
#> 
#> $coef_sign_index$private_cabins_cumulative_[[17]]
#> [1] 5 7
#> 
#> $coef_sign_index$private_cabins_cumulative_[[18]]
#> [1] 5 7
#> 
#> $coef_sign_index$private_cabins_cumulative_[[19]]
#> [1] 5
#> 
#> 
#> $coef_sign_index$public_cabins_high_cumulative_
#> $coef_sign_index$public_cabins_high_cumulative_[[1]]
#> [1] 7
#> 
#> $coef_sign_index$public_cabins_high_cumulative_[[2]]
#> [1] 5 7
#> 
#> $coef_sign_index$public_cabins_high_cumulative_[[3]]
#> [1] 7
#> 
#> $coef_sign_index$public_cabins_high_cumulative_[[4]]
#> [1] 5 7
#> 
#> $coef_sign_index$public_cabins_high_cumulative_[[5]]
#> [1] 7
#> 
#> $coef_sign_index$public_cabins_high_cumulative_[[6]]
#> [1] 7
#> 
#> $coef_sign_index$public_cabins_high_cumulative_[[7]]
#> [1] 7
#> 
#> $coef_sign_index$public_cabins_high_cumulative_[[8]]
#> [1] 3 7
#> 
#> $coef_sign_index$public_cabins_high_cumulative_[[9]]
#> [1] 4 7
#> 
#> $coef_sign_index$public_cabins_high_cumulative_[[10]]
#> [1] 7
#> 
#> $coef_sign_index$public_cabins_high_cumulative_[[11]]
#> [1] 7
#> 
#> $coef_sign_index$public_cabins_high_cumulative_[[12]]
#> [1] 7
#> 
#> $coef_sign_index$public_cabins_high_cumulative_[[13]]
#> [1] 7
#> 
#> $coef_sign_index$public_cabins_high_cumulative_[[14]]
#> [1] 7
#> 
#> $coef_sign_index$public_cabins_high_cumulative_[[15]]
#> [1] 7
#> 
#> $coef_sign_index$public_cabins_high_cumulative_[[16]]
#> [1] 7
#> 
#> $coef_sign_index$public_cabins_high_cumulative_[[17]]
#> [1] 7
#> 
#> $coef_sign_index$public_cabins_high_cumulative_[[18]]
#> [1] 7
#> 
#> $coef_sign_index$public_cabins_high_cumulative_[[19]]
#> [1] 7
#> 
#> 
#> $coef_sign_index$trails_cumulative_
#> $coef_sign_index$trails_cumulative_[[1]]
#> [1] 5 6
#> 
#> $coef_sign_index$trails_cumulative_[[2]]
#> [1] 6
#> 
#> $coef_sign_index$trails_cumulative_[[3]]
#> [1] 6
#> 
#> $coef_sign_index$trails_cumulative_[[4]]
#> [1] 1 6
#> 
#> $coef_sign_index$trails_cumulative_[[5]]
#> [1] 2 3 6
#> 
#> $coef_sign_index$trails_cumulative_[[6]]
#> [1] 6
#> 
#> $coef_sign_index$trails_cumulative_[[7]]
#> [1] 1 6
#> 
#> $coef_sign_index$trails_cumulative_[[8]]
#> [1] 3 6 8
#> 
#> $coef_sign_index$trails_cumulative_[[9]]
#> [1] 6
#> 
#> $coef_sign_index$trails_cumulative_[[10]]
#> [1] 6
#> 
#> $coef_sign_index$trails_cumulative_[[11]]
#> [1] 1 5 6
#> 
#> $coef_sign_index$trails_cumulative_[[12]]
#> [1] 6
#> 
#> $coef_sign_index$trails_cumulative_[[13]]
#> [1] 6
#> 
#> $coef_sign_index$trails_cumulative_[[14]]
#> [1] 1 6
#> 
#> $coef_sign_index$trails_cumulative_[[15]]
#> [1] 2 6
#> 
#> $coef_sign_index$trails_cumulative_[[16]]
#> [1] 6
#> 
#> $coef_sign_index$trails_cumulative_[[17]]
#> [1] 3 6
#> 
#> $coef_sign_index$trails_cumulative_[[18]]
#> [1] 6
#> 
#> $coef_sign_index$trails_cumulative_[[19]]
#> [1] 3 6
#> 
#> 
#> 
#> $coef_sign_names
#> $coef_sign_names$private_cabins_cumulative_
#> $coef_sign_names$private_cabins_cumulative_[[1]]
#> character(0)
#> 
#> $coef_sign_names$private_cabins_cumulative_[[2]]
#> [1] "private_cabins_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$private_cabins_cumulative_[[3]]
#> character(0)
#> 
#> $coef_sign_names$private_cabins_cumulative_[[4]]
#> [1] "private_cabins_cumulative_exp_decay2500" 
#> [2] "private_cabins_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$private_cabins_cumulative_[[5]]
#> [1] "private_cabins_cumulative_exp_decay2500" 
#> [2] "private_cabins_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$private_cabins_cumulative_[[6]]
#> [1] "private_cabins_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$private_cabins_cumulative_[[7]]
#> character(0)
#> 
#> $coef_sign_names$private_cabins_cumulative_[[8]]
#> [1] "private_cabins_cumulative_exp_decay2500" 
#> [2] "private_cabins_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$private_cabins_cumulative_[[9]]
#> [1] "private_cabins_cumulative_exp_decay2500" 
#> [2] "private_cabins_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$private_cabins_cumulative_[[10]]
#> [1] "private_cabins_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$private_cabins_cumulative_[[11]]
#> [1] "private_cabins_cumulative_exp_decay2500" 
#> [2] "private_cabins_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$private_cabins_cumulative_[[12]]
#> character(0)
#> 
#> $coef_sign_names$private_cabins_cumulative_[[13]]
#> [1] "private_cabins_cumulative_exp_decay2500" 
#> [2] "private_cabins_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$private_cabins_cumulative_[[14]]
#> [1] "private_cabins_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$private_cabins_cumulative_[[15]]
#> [1] "private_cabins_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$private_cabins_cumulative_[[16]]
#> [1] "private_cabins_cumulative_exp_decay2500" 
#> [2] "private_cabins_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$private_cabins_cumulative_[[17]]
#> [1] "private_cabins_cumulative_exp_decay2500" 
#> [2] "private_cabins_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$private_cabins_cumulative_[[18]]
#> [1] "private_cabins_cumulative_exp_decay2500" 
#> [2] "private_cabins_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$private_cabins_cumulative_[[19]]
#> [1] "private_cabins_cumulative_exp_decay2500"
#> 
#> 
#> $coef_sign_names$public_cabins_high_cumulative_
#> $coef_sign_names$public_cabins_high_cumulative_[[1]]
#> [1] "public_cabins_high_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$public_cabins_high_cumulative_[[2]]
#> [1] "public_cabins_high_cumulative_exp_decay2500" 
#> [2] "public_cabins_high_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$public_cabins_high_cumulative_[[3]]
#> [1] "public_cabins_high_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$public_cabins_high_cumulative_[[4]]
#> [1] "public_cabins_high_cumulative_exp_decay2500" 
#> [2] "public_cabins_high_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$public_cabins_high_cumulative_[[5]]
#> [1] "public_cabins_high_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$public_cabins_high_cumulative_[[6]]
#> [1] "public_cabins_high_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$public_cabins_high_cumulative_[[7]]
#> [1] "public_cabins_high_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$public_cabins_high_cumulative_[[8]]
#> [1] "public_cabins_high_cumulative_exp_decay500"  
#> [2] "public_cabins_high_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$public_cabins_high_cumulative_[[9]]
#> [1] "public_cabins_high_cumulative_exp_decay1000" 
#> [2] "public_cabins_high_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$public_cabins_high_cumulative_[[10]]
#> [1] "public_cabins_high_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$public_cabins_high_cumulative_[[11]]
#> [1] "public_cabins_high_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$public_cabins_high_cumulative_[[12]]
#> [1] "public_cabins_high_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$public_cabins_high_cumulative_[[13]]
#> [1] "public_cabins_high_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$public_cabins_high_cumulative_[[14]]
#> [1] "public_cabins_high_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$public_cabins_high_cumulative_[[15]]
#> [1] "public_cabins_high_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$public_cabins_high_cumulative_[[16]]
#> [1] "public_cabins_high_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$public_cabins_high_cumulative_[[17]]
#> [1] "public_cabins_high_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$public_cabins_high_cumulative_[[18]]
#> [1] "public_cabins_high_cumulative_exp_decay10000"
#> 
#> $coef_sign_names$public_cabins_high_cumulative_[[19]]
#> [1] "public_cabins_high_cumulative_exp_decay10000"
#> 
#> 
#> $coef_sign_names$trails_cumulative_
#> $coef_sign_names$trails_cumulative_[[1]]
#> [1] "trails_cumulative_exp_decay2500" "trails_cumulative_exp_decay5000"
#> 
#> $coef_sign_names$trails_cumulative_[[2]]
#> [1] "trails_cumulative_exp_decay5000"
#> 
#> $coef_sign_names$trails_cumulative_[[3]]
#> [1] "trails_cumulative_exp_decay5000"
#> 
#> $coef_sign_names$trails_cumulative_[[4]]
#> [1] "trails_cumulative_exp_decay100"  "trails_cumulative_exp_decay5000"
#> 
#> $coef_sign_names$trails_cumulative_[[5]]
#> [1] "trails_cumulative_exp_decay250"  "trails_cumulative_exp_decay500" 
#> [3] "trails_cumulative_exp_decay5000"
#> 
#> $coef_sign_names$trails_cumulative_[[6]]
#> [1] "trails_cumulative_exp_decay5000"
#> 
#> $coef_sign_names$trails_cumulative_[[7]]
#> [1] "trails_cumulative_exp_decay100"  "trails_cumulative_exp_decay5000"
#> 
#> $coef_sign_names$trails_cumulative_[[8]]
#> [1] "trails_cumulative_exp_decay500"   "trails_cumulative_exp_decay5000" 
#> [3] "trails_cumulative_exp_decay20000"
#> 
#> $coef_sign_names$trails_cumulative_[[9]]
#> [1] "trails_cumulative_exp_decay5000"
#> 
#> $coef_sign_names$trails_cumulative_[[10]]
#> [1] "trails_cumulative_exp_decay5000"
#> 
#> $coef_sign_names$trails_cumulative_[[11]]
#> [1] "trails_cumulative_exp_decay100"  "trails_cumulative_exp_decay2500"
#> [3] "trails_cumulative_exp_decay5000"
#> 
#> $coef_sign_names$trails_cumulative_[[12]]
#> [1] "trails_cumulative_exp_decay5000"
#> 
#> $coef_sign_names$trails_cumulative_[[13]]
#> [1] "trails_cumulative_exp_decay5000"
#> 
#> $coef_sign_names$trails_cumulative_[[14]]
#> [1] "trails_cumulative_exp_decay100"  "trails_cumulative_exp_decay5000"
#> 
#> $coef_sign_names$trails_cumulative_[[15]]
#> [1] "trails_cumulative_exp_decay250"  "trails_cumulative_exp_decay5000"
#> 
#> $coef_sign_names$trails_cumulative_[[16]]
#> [1] "trails_cumulative_exp_decay5000"
#> 
#> $coef_sign_names$trails_cumulative_[[17]]
#> [1] "trails_cumulative_exp_decay500"  "trails_cumulative_exp_decay5000"
#> 
#> $coef_sign_names$trails_cumulative_[[18]]
#> [1] "trails_cumulative_exp_decay5000"
#> 
#> $coef_sign_names$trails_cumulative_[[19]]
#> [1] "trails_cumulative_exp_decay500"  "trails_cumulative_exp_decay5000"
#> 
#> 
#> 
#> $coef_sign_radii
#> $coef_sign_radii$private_cabins_cumulative_
#> $coef_sign_radii$private_cabins_cumulative_[[1]]
#> numeric(0)
#> 
#> $coef_sign_radii$private_cabins_cumulative_[[2]]
#> [1] 10000
#> 
#> $coef_sign_radii$private_cabins_cumulative_[[3]]
#> numeric(0)
#> 
#> $coef_sign_radii$private_cabins_cumulative_[[4]]
#> [1]  2500 10000
#> 
#> $coef_sign_radii$private_cabins_cumulative_[[5]]
#> [1]  2500 10000
#> 
#> $coef_sign_radii$private_cabins_cumulative_[[6]]
#> [1] 10000
#> 
#> $coef_sign_radii$private_cabins_cumulative_[[7]]
#> numeric(0)
#> 
#> $coef_sign_radii$private_cabins_cumulative_[[8]]
#> [1]  2500 10000
#> 
#> $coef_sign_radii$private_cabins_cumulative_[[9]]
#> [1]  2500 10000
#> 
#> $coef_sign_radii$private_cabins_cumulative_[[10]]
#> [1] 10000
#> 
#> $coef_sign_radii$private_cabins_cumulative_[[11]]
#> [1]  2500 10000
#> 
#> $coef_sign_radii$private_cabins_cumulative_[[12]]
#> numeric(0)
#> 
#> $coef_sign_radii$private_cabins_cumulative_[[13]]
#> [1]  2500 10000
#> 
#> $coef_sign_radii$private_cabins_cumulative_[[14]]
#> [1] 10000
#> 
#> $coef_sign_radii$private_cabins_cumulative_[[15]]
#> [1] 10000
#> 
#> $coef_sign_radii$private_cabins_cumulative_[[16]]
#> [1]  2500 10000
#> 
#> $coef_sign_radii$private_cabins_cumulative_[[17]]
#> [1]  2500 10000
#> 
#> $coef_sign_radii$private_cabins_cumulative_[[18]]
#> [1]  2500 10000
#> 
#> $coef_sign_radii$private_cabins_cumulative_[[19]]
#> [1] 2500
#> 
#> 
#> $coef_sign_radii$public_cabins_high_cumulative_
#> $coef_sign_radii$public_cabins_high_cumulative_[[1]]
#> [1] 10000
#> 
#> $coef_sign_radii$public_cabins_high_cumulative_[[2]]
#> [1]  2500 10000
#> 
#> $coef_sign_radii$public_cabins_high_cumulative_[[3]]
#> [1] 10000
#> 
#> $coef_sign_radii$public_cabins_high_cumulative_[[4]]
#> [1]  2500 10000
#> 
#> $coef_sign_radii$public_cabins_high_cumulative_[[5]]
#> [1] 10000
#> 
#> $coef_sign_radii$public_cabins_high_cumulative_[[6]]
#> [1] 10000
#> 
#> $coef_sign_radii$public_cabins_high_cumulative_[[7]]
#> [1] 10000
#> 
#> $coef_sign_radii$public_cabins_high_cumulative_[[8]]
#> [1]   500 10000
#> 
#> $coef_sign_radii$public_cabins_high_cumulative_[[9]]
#> [1]  1000 10000
#> 
#> $coef_sign_radii$public_cabins_high_cumulative_[[10]]
#> [1] 10000
#> 
#> $coef_sign_radii$public_cabins_high_cumulative_[[11]]
#> [1] 10000
#> 
#> $coef_sign_radii$public_cabins_high_cumulative_[[12]]
#> [1] 10000
#> 
#> $coef_sign_radii$public_cabins_high_cumulative_[[13]]
#> [1] 10000
#> 
#> $coef_sign_radii$public_cabins_high_cumulative_[[14]]
#> [1] 10000
#> 
#> $coef_sign_radii$public_cabins_high_cumulative_[[15]]
#> [1] 10000
#> 
#> $coef_sign_radii$public_cabins_high_cumulative_[[16]]
#> [1] 10000
#> 
#> $coef_sign_radii$public_cabins_high_cumulative_[[17]]
#> [1] 10000
#> 
#> $coef_sign_radii$public_cabins_high_cumulative_[[18]]
#> [1] 10000
#> 
#> $coef_sign_radii$public_cabins_high_cumulative_[[19]]
#> [1] 10000
#> 
#> 
#> $coef_sign_radii$trails_cumulative_
#> $coef_sign_radii$trails_cumulative_[[1]]
#> [1] 2500 5000
#> 
#> $coef_sign_radii$trails_cumulative_[[2]]
#> [1] 5000
#> 
#> $coef_sign_radii$trails_cumulative_[[3]]
#> [1] 5000
#> 
#> $coef_sign_radii$trails_cumulative_[[4]]
#> [1]  100 5000
#> 
#> $coef_sign_radii$trails_cumulative_[[5]]
#> [1]  250  500 5000
#> 
#> $coef_sign_radii$trails_cumulative_[[6]]
#> [1] 5000
#> 
#> $coef_sign_radii$trails_cumulative_[[7]]
#> [1]  100 5000
#> 
#> $coef_sign_radii$trails_cumulative_[[8]]
#> [1]   500  5000 20000
#> 
#> $coef_sign_radii$trails_cumulative_[[9]]
#> [1] 5000
#> 
#> $coef_sign_radii$trails_cumulative_[[10]]
#> [1] 5000
#> 
#> $coef_sign_radii$trails_cumulative_[[11]]
#> [1]  100 2500 5000
#> 
#> $coef_sign_radii$trails_cumulative_[[12]]
#> [1] 5000
#> 
#> $coef_sign_radii$trails_cumulative_[[13]]
#> [1] 5000
#> 
#> $coef_sign_radii$trails_cumulative_[[14]]
#> [1]  100 5000
#> 
#> $coef_sign_radii$trails_cumulative_[[15]]
#> [1]  250 5000
#> 
#> $coef_sign_radii$trails_cumulative_[[16]]
#> [1] 5000
#> 
#> $coef_sign_radii$trails_cumulative_[[17]]
#> [1]  500 5000
#> 
#> $coef_sign_radii$trails_cumulative_[[18]]
#> [1] 5000
#> 
#> $coef_sign_radii$trails_cumulative_[[19]]
#> [1]  500 5000
#> 
#> 
#> 
#> $coef_sign_value
#> $coef_sign_value$private_cabins_cumulative_
#> $coef_sign_value$private_cabins_cumulative_[[1]]
#> numeric(0)
#> 
#> $coef_sign_value$private_cabins_cumulative_[[2]]
#> [1] -0.01493428
#> 
#> $coef_sign_value$private_cabins_cumulative_[[3]]
#> numeric(0)
#> 
#> $coef_sign_value$private_cabins_cumulative_[[4]]
#> [1]  0.00000000 -0.01493428
#> 
#> $coef_sign_value$private_cabins_cumulative_[[5]]
#> [1]  0.00000000 -0.01493428
#> 
#> $coef_sign_value$private_cabins_cumulative_[[6]]
#> [1] -0.01493428
#> 
#> $coef_sign_value$private_cabins_cumulative_[[7]]
#> numeric(0)
#> 
#> $coef_sign_value$private_cabins_cumulative_[[8]]
#> [1]  0.00000000 -0.01493428
#> 
#> $coef_sign_value$private_cabins_cumulative_[[9]]
#> [1]  0.00000000 -0.01493428
#> 
#> $coef_sign_value$private_cabins_cumulative_[[10]]
#> [1] -0.01493428
#> 
#> $coef_sign_value$private_cabins_cumulative_[[11]]
#> [1]  0.00000000 -0.01493428
#> 
#> $coef_sign_value$private_cabins_cumulative_[[12]]
#> numeric(0)
#> 
#> $coef_sign_value$private_cabins_cumulative_[[13]]
#> [1]  0.00000000 -0.01493428
#> 
#> $coef_sign_value$private_cabins_cumulative_[[14]]
#> [1] -0.01493428
#> 
#> $coef_sign_value$private_cabins_cumulative_[[15]]
#> [1] -0.01493428
#> 
#> $coef_sign_value$private_cabins_cumulative_[[16]]
#> [1]  0.00000000 -0.01493428
#> 
#> $coef_sign_value$private_cabins_cumulative_[[17]]
#> [1]  0.00000000 -0.01493428
#> 
#> $coef_sign_value$private_cabins_cumulative_[[18]]
#> [1]  0.00000000 -0.01493428
#> 
#> $coef_sign_value$private_cabins_cumulative_[[19]]
#> [1] 0
#> 
#> 
#> $coef_sign_value$public_cabins_high_cumulative_
#> $coef_sign_value$public_cabins_high_cumulative_[[1]]
#> [1] 1.584815
#> 
#> $coef_sign_value$public_cabins_high_cumulative_[[2]]
#> [1] 0.000000 1.584815
#> 
#> $coef_sign_value$public_cabins_high_cumulative_[[3]]
#> [1] 1.584815
#> 
#> $coef_sign_value$public_cabins_high_cumulative_[[4]]
#> [1] 0.000000 1.584815
#> 
#> $coef_sign_value$public_cabins_high_cumulative_[[5]]
#> [1] 1.584815
#> 
#> $coef_sign_value$public_cabins_high_cumulative_[[6]]
#> [1] 1.584815
#> 
#> $coef_sign_value$public_cabins_high_cumulative_[[7]]
#> [1] 1.584815
#> 
#> $coef_sign_value$public_cabins_high_cumulative_[[8]]
#> [1] 0.000000 1.584815
#> 
#> $coef_sign_value$public_cabins_high_cumulative_[[9]]
#> [1] 0.000000 1.584815
#> 
#> $coef_sign_value$public_cabins_high_cumulative_[[10]]
#> [1] 1.584815
#> 
#> $coef_sign_value$public_cabins_high_cumulative_[[11]]
#> [1] 1.584815
#> 
#> $coef_sign_value$public_cabins_high_cumulative_[[12]]
#> [1] 1.584815
#> 
#> $coef_sign_value$public_cabins_high_cumulative_[[13]]
#> [1] 1.584815
#> 
#> $coef_sign_value$public_cabins_high_cumulative_[[14]]
#> [1] 1.584815
#> 
#> $coef_sign_value$public_cabins_high_cumulative_[[15]]
#> [1] 1.584815
#> 
#> $coef_sign_value$public_cabins_high_cumulative_[[16]]
#> [1] 1.584815
#> 
#> $coef_sign_value$public_cabins_high_cumulative_[[17]]
#> [1] 1.584815
#> 
#> $coef_sign_value$public_cabins_high_cumulative_[[18]]
#> [1] 1.584815
#> 
#> $coef_sign_value$public_cabins_high_cumulative_[[19]]
#> [1] 1.584815
#> 
#> 
#> $coef_sign_value$trails_cumulative_
#> $coef_sign_value$trails_cumulative_[[1]]
#> [1] 0.02547249 0.01438307
#> 
#> $coef_sign_value$trails_cumulative_[[2]]
#> [1] 0.01438307
#> 
#> $coef_sign_value$trails_cumulative_[[3]]
#> [1] 0.01438307
#> 
#> $coef_sign_value$trails_cumulative_[[4]]
#> [1] 0.00000000 0.01438307
#> 
#> $coef_sign_value$trails_cumulative_[[5]]
#> [1]  0.00000000 -0.24583477  0.01438307
#> 
#> $coef_sign_value$trails_cumulative_[[6]]
#> [1] 0.01438307
#> 
#> $coef_sign_value$trails_cumulative_[[7]]
#> [1] 0.00000000 0.01438307
#> 
#> $coef_sign_value$trails_cumulative_[[8]]
#> [1] -0.245834772  0.014383067 -0.005513045
#> 
#> $coef_sign_value$trails_cumulative_[[9]]
#> [1] 0.01438307
#> 
#> $coef_sign_value$trails_cumulative_[[10]]
#> [1] 0.01438307
#> 
#> $coef_sign_value$trails_cumulative_[[11]]
#> [1] 0.00000000 0.02547249 0.01438307
#> 
#> $coef_sign_value$trails_cumulative_[[12]]
#> [1] 0.01438307
#> 
#> $coef_sign_value$trails_cumulative_[[13]]
#> [1] 0.01438307
#> 
#> $coef_sign_value$trails_cumulative_[[14]]
#> [1] 0.00000000 0.01438307
#> 
#> $coef_sign_value$trails_cumulative_[[15]]
#> [1] 0.00000000 0.01438307
#> 
#> $coef_sign_value$trails_cumulative_[[16]]
#> [1] 0.01438307
#> 
#> $coef_sign_value$trails_cumulative_[[17]]
#> [1] -0.24583477  0.01438307
#> 
#> $coef_sign_value$trails_cumulative_[[18]]
#> [1] 0.01438307
#> 
#> $coef_sign_value$trails_cumulative_[[19]]
#> [1] -0.24583477  0.01438307
#> 
#> 
#> 
#> $coef_sign
#>       private_cabins_cumulative_ public_cabins_high_cumulative_
#>  [1,]                          0                              1
#>  [2,]                          1                              2
#>  [3,]                          0                              1
#>  [4,]                          2                              2
#>  [5,]                          2                              1
#>  [6,]                          1                              1
#>  [7,]                          0                              1
#>  [8,]                          2                              2
#>  [9,]                          2                              2
#> [10,]                          1                              1
#> [11,]                          2                              1
#> [12,]                          0                              1
#> [13,]                          2                              1
#> [14,]                          1                              1
#> [15,]                          1                              1
#> [16,]                          2                              1
#> [17,]                          2                              1
#> [18,]                          2                              1
#> [19,]                          1                              1
#>       trails_cumulative_
#>  [1,]                  2
#>  [2,]                  1
#>  [3,]                  1
#>  [4,]                  2
#>  [5,]                  3
#>  [6,]                  1
#>  [7,]                  2
#>  [8,]                  3
#>  [9,]                  1
#> [10,]                  1
#> [11,]                  3
#> [12,]                  1
#> [13,]                  1
#> [14,]                  2
#> [15,]                  2
#> [16,]                  1
#> [17,]                  2
#> [18,]                  1
#> [19,]                  2
#> 
#> $coef_sign_sum
#> [1] 79
#> 
#> $cross_index
#> $cross_index$private_cabins_cumulative_
#> $cross_index$private_cabins_cumulative_$Resample01
#> integer(0)
#> 
#> $cross_index$private_cabins_cumulative_$Resample02
#> integer(0)
#> 
#> $cross_index$private_cabins_cumulative_$Resample03
#> integer(0)
#> 
#> $cross_index$private_cabins_cumulative_$Resample04
#> integer(0)
#> 
#> $cross_index$private_cabins_cumulative_$Resample05
#> integer(0)
#> 
#> $cross_index$private_cabins_cumulative_$Resample06
#> integer(0)
#> 
#> $cross_index$private_cabins_cumulative_$Resample08
#> integer(0)
#> 
#> $cross_index$private_cabins_cumulative_$Resample09
#> integer(0)
#> 
#> $cross_index$private_cabins_cumulative_$Resample10
#> integer(0)
#> 
#> $cross_index$private_cabins_cumulative_$Resample11
#> integer(0)
#> 
#> $cross_index$private_cabins_cumulative_$Resample12
#> integer(0)
#> 
#> $cross_index$private_cabins_cumulative_$Resample13
#> integer(0)
#> 
#> $cross_index$private_cabins_cumulative_$Resample14
#> integer(0)
#> 
#> $cross_index$private_cabins_cumulative_$Resample15
#> integer(0)
#> 
#> $cross_index$private_cabins_cumulative_$Resample16
#> integer(0)
#> 
#> $cross_index$private_cabins_cumulative_$Resample17
#> integer(0)
#> 
#> $cross_index$private_cabins_cumulative_$Resample18
#> integer(0)
#> 
#> $cross_index$private_cabins_cumulative_$Resample19
#> integer(0)
#> 
#> $cross_index$private_cabins_cumulative_$Resample20
#> integer(0)
#> 
#> 
#> $cross_index$public_cabins_high_cumulative_
#> $cross_index$public_cabins_high_cumulative_$Resample01
#> integer(0)
#> 
#> $cross_index$public_cabins_high_cumulative_$Resample02
#> integer(0)
#> 
#> $cross_index$public_cabins_high_cumulative_$Resample03
#> integer(0)
#> 
#> $cross_index$public_cabins_high_cumulative_$Resample04
#> integer(0)
#> 
#> $cross_index$public_cabins_high_cumulative_$Resample05
#> integer(0)
#> 
#> $cross_index$public_cabins_high_cumulative_$Resample06
#> integer(0)
#> 
#> $cross_index$public_cabins_high_cumulative_$Resample08
#> integer(0)
#> 
#> $cross_index$public_cabins_high_cumulative_$Resample09
#> [1]  46 182
#> 
#> $cross_index$public_cabins_high_cumulative_$Resample10
#> [1]  53 248
#> 
#> $cross_index$public_cabins_high_cumulative_$Resample11
#> integer(0)
#> 
#> $cross_index$public_cabins_high_cumulative_$Resample12
#> integer(0)
#> 
#> $cross_index$public_cabins_high_cumulative_$Resample13
#> integer(0)
#> 
#> $cross_index$public_cabins_high_cumulative_$Resample14
#> integer(0)
#> 
#> $cross_index$public_cabins_high_cumulative_$Resample15
#> integer(0)
#> 
#> $cross_index$public_cabins_high_cumulative_$Resample16
#> integer(0)
#> 
#> $cross_index$public_cabins_high_cumulative_$Resample17
#> integer(0)
#> 
#> $cross_index$public_cabins_high_cumulative_$Resample18
#> integer(0)
#> 
#> $cross_index$public_cabins_high_cumulative_$Resample19
#> integer(0)
#> 
#> $cross_index$public_cabins_high_cumulative_$Resample20
#> integer(0)
#> 
#> 
#> $cross_index$trails_cumulative_
#> $cross_index$trails_cumulative_$Resample01
#> [1]  203 1389
#> 
#> $cross_index$trails_cumulative_$Resample02
#> [1]  294 1777
#> 
#> $cross_index$trails_cumulative_$Resample03
#> [1]   59 1496
#> 
#> $cross_index$trails_cumulative_$Resample04
#> [1]   67  172 1605
#> 
#> $cross_index$trails_cumulative_$Resample05
#> [1]  106  401 2011
#> 
#> $cross_index$trails_cumulative_$Resample06
#> [1]  266 1658
#> 
#> $cross_index$trails_cumulative_$Resample08
#> [1]    7  238 1634
#> 
#> $cross_index$trails_cumulative_$Resample09
#> [1]  397 1603
#> 
#> $cross_index$trails_cumulative_$Resample10
#> [1]  327 1899
#> 
#> $cross_index$trails_cumulative_$Resample11
#> [1]  356 1822
#> 
#> $cross_index$trails_cumulative_$Resample12
#> [1]  261 1560
#> 
#> $cross_index$trails_cumulative_$Resample13
#> [1]  236 1932
#> 
#> $cross_index$trails_cumulative_$Resample14
#> [1]  268 1709
#> 
#> $cross_index$trails_cumulative_$Resample15
#> [1]   11  283 1914
#> 
#> $cross_index$trails_cumulative_$Resample16
#> [1]   33  267 1679
#> 
#> $cross_index$trails_cumulative_$Resample17
#> [1]  454 1906
#> 
#> $cross_index$trails_cumulative_$Resample18
#> [1]  387 1794
#> 
#> $cross_index$trails_cumulative_$Resample19
#> [1]  335 1566
#> 
#> $cross_index$trails_cumulative_$Resample20
#> [1]  379 1877
#> 
#> 
#> 
#> $where_crosses
#> $where_crosses$private_cabins_cumulative_
#> $where_crosses$private_cabins_cumulative_[[1]]
#> numeric(0)
#> 
#> $where_crosses$private_cabins_cumulative_[[2]]
#> numeric(0)
#> 
#> $where_crosses$private_cabins_cumulative_[[3]]
#> numeric(0)
#> 
#> $where_crosses$private_cabins_cumulative_[[4]]
#> numeric(0)
#> 
#> $where_crosses$private_cabins_cumulative_[[5]]
#> numeric(0)
#> 
#> $where_crosses$private_cabins_cumulative_[[6]]
#> numeric(0)
#> 
#> $where_crosses$private_cabins_cumulative_[[7]]
#> numeric(0)
#> 
#> $where_crosses$private_cabins_cumulative_[[8]]
#> numeric(0)
#> 
#> $where_crosses$private_cabins_cumulative_[[9]]
#> numeric(0)
#> 
#> $where_crosses$private_cabins_cumulative_[[10]]
#> numeric(0)
#> 
#> $where_crosses$private_cabins_cumulative_[[11]]
#> numeric(0)
#> 
#> $where_crosses$private_cabins_cumulative_[[12]]
#> numeric(0)
#> 
#> $where_crosses$private_cabins_cumulative_[[13]]
#> numeric(0)
#> 
#> $where_crosses$private_cabins_cumulative_[[14]]
#> numeric(0)
#> 
#> $where_crosses$private_cabins_cumulative_[[15]]
#> numeric(0)
#> 
#> $where_crosses$private_cabins_cumulative_[[16]]
#> numeric(0)
#> 
#> $where_crosses$private_cabins_cumulative_[[17]]
#> numeric(0)
#> 
#> $where_crosses$private_cabins_cumulative_[[18]]
#> numeric(0)
#> 
#> $where_crosses$private_cabins_cumulative_[[19]]
#> numeric(0)
#> 
#> 
#> $where_crosses$public_cabins_high_cumulative_
#> $where_crosses$public_cabins_high_cumulative_[[1]]
#> numeric(0)
#> 
#> $where_crosses$public_cabins_high_cumulative_[[2]]
#> numeric(0)
#> 
#> $where_crosses$public_cabins_high_cumulative_[[3]]
#> numeric(0)
#> 
#> $where_crosses$public_cabins_high_cumulative_[[4]]
#> numeric(0)
#> 
#> $where_crosses$public_cabins_high_cumulative_[[5]]
#> numeric(0)
#> 
#> $where_crosses$public_cabins_high_cumulative_[[6]]
#> numeric(0)
#> 
#> $where_crosses$public_cabins_high_cumulative_[[7]]
#> numeric(0)
#> 
#> $where_crosses$public_cabins_high_cumulative_[[8]]
#> [1]  90 362
#> 
#> $where_crosses$public_cabins_high_cumulative_[[9]]
#> [1] 104 494
#> 
#> $where_crosses$public_cabins_high_cumulative_[[10]]
#> numeric(0)
#> 
#> $where_crosses$public_cabins_high_cumulative_[[11]]
#> numeric(0)
#> 
#> $where_crosses$public_cabins_high_cumulative_[[12]]
#> numeric(0)
#> 
#> $where_crosses$public_cabins_high_cumulative_[[13]]
#> numeric(0)
#> 
#> $where_crosses$public_cabins_high_cumulative_[[14]]
#> numeric(0)
#> 
#> $where_crosses$public_cabins_high_cumulative_[[15]]
#> numeric(0)
#> 
#> $where_crosses$public_cabins_high_cumulative_[[16]]
#> numeric(0)
#> 
#> $where_crosses$public_cabins_high_cumulative_[[17]]
#> numeric(0)
#> 
#> $where_crosses$public_cabins_high_cumulative_[[18]]
#> numeric(0)
#> 
#> $where_crosses$public_cabins_high_cumulative_[[19]]
#> numeric(0)
#> 
#> 
#> $where_crosses$trails_cumulative_
#> $where_crosses$trails_cumulative_[[1]]
#> [1]  404 2776
#> 
#> $where_crosses$trails_cumulative_[[2]]
#> [1]  586 3552
#> 
#> $where_crosses$trails_cumulative_[[3]]
#> [1]  116 2990
#> 
#> $where_crosses$trails_cumulative_[[4]]
#> [1]  132  342 3208
#> 
#> $where_crosses$trails_cumulative_[[5]]
#> [1]  210  800 4020
#> 
#> $where_crosses$trails_cumulative_[[6]]
#> [1]  530 3314
#> 
#> $where_crosses$trails_cumulative_[[7]]
#> [1]   12  474 3266
#> 
#> $where_crosses$trails_cumulative_[[8]]
#> [1]  792 3204
#> 
#> $where_crosses$trails_cumulative_[[9]]
#> [1]  652 3796
#> 
#> $where_crosses$trails_cumulative_[[10]]
#> [1]  710 3642
#> 
#> $where_crosses$trails_cumulative_[[11]]
#> [1]  520 3118
#> 
#> $where_crosses$trails_cumulative_[[12]]
#> [1]  470 3862
#> 
#> $where_crosses$trails_cumulative_[[13]]
#> [1]  534 3416
#> 
#> $where_crosses$trails_cumulative_[[14]]
#> [1]   20  564 3826
#> 
#> $where_crosses$trails_cumulative_[[15]]
#> [1]   64  532 3356
#> 
#> $where_crosses$trails_cumulative_[[16]]
#> [1]  906 3810
#> 
#> $where_crosses$trails_cumulative_[[17]]
#> [1]  772 3586
#> 
#> $where_crosses$trails_cumulative_[[18]]
#> [1]  668 3130
#> 
#> $where_crosses$trails_cumulative_[[19]]
#> [1]  756 3752
#> 
#> 
#> 
#> $n_crosses
#>            private_cabins_cumulative_ public_cabins_high_cumulative_
#> Resample01                          0                              0
#> Resample02                          0                              0
#> Resample03                          0                              0
#> Resample04                          0                              0
#> Resample05                          0                              0
#> Resample06                          0                              0
#> Resample08                          0                              0
#> Resample09                          0                              2
#> Resample10                          0                              2
#> Resample11                          0                              0
#> Resample12                          0                              0
#> Resample13                          0                              0
#> Resample14                          0                              0
#> Resample15                          0                              0
#> Resample16                          0                              0
#> Resample17                          0                              0
#> Resample18                          0                              0
#> Resample19                          0                              0
#> Resample20                          0                              0
#>            trails_cumulative_
#> Resample01                  2
#> Resample02                  2
#> Resample03                  2
#> Resample04                  3
#> Resample05                  3
#> Resample06                  2
#> Resample08                  3
#> Resample09                  2
#> Resample10                  2
#> Resample11                  2
#> Resample12                  2
#> Resample13                  2
#> Resample14                  2
#> Resample15                  3
#> Resample16                  3
#> Resample17                  2
#> Resample18                  2
#> Resample19                  2
#> Resample20                  2
#> 
#> $n_crosses_total
#> [1] 47
#> 
#> $response_area_opposite
#> [1] "This needs to be implemented for individual models. Please raise na issue on our Github repo."
#> 
#> $response_area_opposite_total
#> [1] "This needs to be implemented for individual models. Please raise na issue on our Github repo."
#> 
#> $response_area_ratio
#> [1] "This needs to be implemented for individual models. Please raise na issue on our Github repo."
#> 
#> $response_area_ratio_total
#> [1] "This needs to be implemented for individual models. Please raise na issue on our Github repo."
#> 
#> $n_inflection
#>     private_cabins_cumulative_ public_cabins_high_cumulative_ 
#>                      0.0000000                      0.2105263 
#>             trails_cumulative_ 
#>                      2.2631579 
#> 
#> $n_inflection_total
#> [1] 2.473684
#> 
#> $difference_inflection
#>     private_cabins_cumulative_ public_cabins_high_cumulative_ 
#>                    0.001461001                    0.462528046 
#>             trails_cumulative_ 
#>                    0.034489552 
#> 
#> $difference_inflection_total
#> [1] 0.4984786
#> 
```
