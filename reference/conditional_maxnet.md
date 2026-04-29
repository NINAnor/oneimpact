# Fits a conditional maxnet model for presence-only species distribution modeling

This function fits a conditional maxnet model using glmnet with
stratified sampling, suitable for presence-only species distribution
models with stratified data structure.

## Usage

``` r
conditional_maxnet(
  p,
  strata,
  data,
  f = maxnet::maxnet.formula(p, data),
  regmult = 1,
  regfun = maxnet::maxnet.default.regularization,
  ...
)
```

## Arguments

- p:

  `[numeric]`  
  Presence indicator vector (1 for presence, 0 for background).

- strata:

  `[factor or numeric]`  
  Stratification factor for conditional logistic regression.

- data:

  `[data.frame]`  
  Data frame containing predictor variables and strata information.

- f:

  `[formula]`  
  Model formula. If not provided, uses
  `maxnet::maxnet.formula(p, data)`.

- regmult:

  `[numeric(1)=1.0]`  
  Regularization multiplier for penalty factors.

- regfun:

  `[function]`  
  Regularization function. Default is
  [`maxnet::maxnet.default.regularization`](https://rdrr.io/pkg/maxnet/man/maxnet.html).

- ...:

  Additional arguments passed to
  [`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html).

## Value

A `maxnet` object with fitted model coefficients, penalties, and data
characteristics.

## See also

[`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html),
[`maxnet::maxnet()`](https://rdrr.io/pkg/maxnet/man/maxnet.html)
