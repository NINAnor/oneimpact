# Prediction based only on step length terms

Computes model predictions based exclusively on movement kernel
variables (e.g. step length and turning angle), using the fitted
coefficients from a conditional logistic regression model.

## Usage

``` r
kernel_prediction(f, data, kernel_vars = c("step_length", "ta"), coefs)
```

## Arguments

- f:

  `[formula]`  
  Model formula, used to identify kernel variable terms.

- data:

  `[data.frame]`  
  Data set on which predictions are made.

- kernel_vars:

  `[character=c("step_length","ta")]`  
  Names of the movement kernel variables in the model.

- coefs:

  `[numeric]`  
  Named numeric vector of fitted model coefficients.

## Value

A numeric matrix of predicted values based on kernel variables only.

## See also

[`fit_net_clogit()`](https://ninanor.github.io/oneimpact/reference/fit_net_functions.md),
[`bag_fit_net_clogit()`](https://ninanor.github.io/oneimpact/reference/bag_fit_net_functions.md)
