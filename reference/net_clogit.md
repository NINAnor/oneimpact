# Fits a conditional logistic regression/SSF/iSSF using glmnet

Low-level wrapper that sets up the design matrix and stratified survival
response for penalized conditional logistic regression, then calls
[`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html)
or
[`glmnet::cv.glmnet()`](https://glmnet.stanford.edu/reference/cv.glmnet.html)
with `family = "cox"`. This function is used internally by
[`fit_net_clogit()`](https://ninanor.github.io/oneimpact/reference/fit_net_functions.md)
and is typically not called directly by the user.

## Usage

``` r
net_clogit(
  f,
  data,
  alpha = 1,
  penalty.factor = NULL,
  type.measure = "deviance",
  standardize = TRUE,
  na.action = "na.pass",
  func = c("glmnet", "cv.glmnet")[1],
  ...
)

net_ssf(
  f,
  data,
  alpha = 1,
  penalty.factor = NULL,
  type.measure = "deviance",
  standardize = TRUE,
  na.action = "na.pass",
  func = c("glmnet", "cv.glmnet")[1],
  ...
)

net_issf(
  f,
  data,
  alpha = 1,
  penalty.factor = NULL,
  type.measure = "deviance",
  standardize = TRUE,
  na.action = "na.pass",
  func = c("glmnet", "cv.glmnet")[1],
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

- alpha:

  `[numeric(1)=1]`  
  Mixing parameter for glmnet. Default is L1-regularization (Lasso),
  with `alpha = 1`. L2-regularization (Ridge regression) is done with
  `alpha = 0`, and elastic-net regression is performed for any `alpha`
  value between `0` and `1`. For more details, see the
  [`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html)
  documentation. For Adaptive and Decay Adaptive Lasso, keep
  `alpha = 1`.

- penalty.factor:

  `[numeric,vector=NULL]`  
  Vector of penalty factors to be used for Adaptive Lasso fitting. The
  vector might have the same length as the the number of columns given
  by the model matrix, `model.matrix(f, data)`. Default is `NULL`, in
  case the same penalty is applied to all variables.

- type.measure:

  `[character(1)="deviance"]`  
  Type of measure to evaluate the model internally in
  [`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html).
  For logistic and conditional logistic regression, it is by default
  `"deviance"`.

- standardize:

  `[logical(1)=TRUE]`  
  Whether glmnet should internally standardize variables before fitting.
  Default is `TRUE`. Set to `FALSE` if variables are already
  standardized.

- na.action:

  `[character(1)="na.pass"]`  
  Default is `"na.pass"`, i.e. rows with NAs are not automatically
  removed from the `model.matrix` used for fitting.

- func:

  `[character(1)="glmnet"]{"glmnet", "cv.glmnet"}`  
  The function to be used for fitting. Default is
  [`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html).
  The second option is
  [`glmnet::cv.glmnet()`](https://glmnet.stanford.edu/reference/cv.glmnet.html)
  which already performs the cross-validation and might include the
  variable selection/callibration within.

- ...:

  `[any]`  
  Additional arguments passed to
  [`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html)
  or
  [`glmnet::cv.glmnet()`](https://glmnet.stanford.edu/reference/cv.glmnet.html).
  Note the `parallel = TRUE` option from glmnet can be passed here.

## Value

A fitted
[`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html)
or
[`glmnet::cv.glmnet()`](https://glmnet.stanford.edu/reference/cv.glmnet.html)
object.

## See also

[`fit_net_clogit()`](https://ninanor.github.io/oneimpact/reference/fit_net_functions.md),
[`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html)
