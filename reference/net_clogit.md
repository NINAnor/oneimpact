# Fits a conditional logistic regression/SSF/iSSF using glmnet

Fits a conditional logistic regression/SSF/iSSF using glmnet

Function with similar name

Function with similar name

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

  Default is L1-regularization (Lasso regression), with `alpha = 1`.
  L2-regularization (Ridge regression) is done with `alpha = 0`, and
  elastic-net regression is performed for any `alpha` value between `0`
  and `1`. For more details, see the
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

  Check option parallel = TRUE from glmnet.
