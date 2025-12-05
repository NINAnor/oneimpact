# Rescale standardized coefficients back to their original range after model fitting

Predictor variables are often standardized to be included in statistical
models and allow comparison of the effect sizes for different
predictors. This functions scales the fitted models coefficients back to
the original scale of the predictors, to allow ecological
interpretation.

## Usage

``` r
rescale_coefficients(...)

# S3 method for class 'coxph'
rescale_coefficients(model, data, ...)

# S3 method for class 'lm'
rescale_coefficients(model, data, ...)

# S3 method for class 'glm'
rescale_coefficients(model, data, ...)

# S3 method for class 'bag'
rescale_coefficients(bag, data, tostd = TRUE, ...)
```

## Arguments

- model:

    
  Fitted model, i.e. the object created by a model fit funtion such as
  "lm", "glm", or "coxph".

- data:

    
  The original `data.frame` with the data used to fit the model.

- bag:

    
  A bag of models, as result of the
  [`bag_models()`](https://ninanor.github.io/oneimpact/reference/bag_models.md)
  function.

- standardize:

  `[logical(1)=TRUE]` If `TRUE` (Default), the coefficients are
  standardized. If `FALSE`, the coefficients are standardized. Only
  numeric coefficients are standardized.

## Value

A vector of rescaled coefficients for the input model.

## Examples

``` r
library(dplyr)

# standardize predictors
iris_std <- iris |>
  dplyr::mutate(across(2:4, ~ scale(.x)))
# fit model
m1 <- lm(Sepal.Length ~ Petal.Length + Species, data = iris_std)
summary(m1)
#> 
#> Call:
#> lm(formula = Sepal.Length ~ Petal.Length + Species, data = iris_std)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -0.75310 -0.23142 -0.00081  0.23085  1.03100 
#> 
#> Coefficients:
#>                   Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)         7.0829     0.1562  45.333  < 2e-16 ***
#> Petal.Length        1.5968     0.1144  13.962  < 2e-16 ***
#> Speciesversicolor  -1.6010     0.1935  -8.275 7.37e-14 ***
#> Speciesvirginica   -2.1177     0.2735  -7.744 1.48e-12 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual standard error: 0.338 on 146 degrees of freedom
#> Multiple R-squared:  0.8367, Adjusted R-squared:  0.8334 
#> F-statistic: 249.4 on 3 and 146 DF,  p-value: < 2.2e-16
#> 

# rescale coefficients
(resc_cf <- rescale_coefficients(m1, iris))
#>       (Intercept)      Petal.Length Speciesversicolor  Speciesvirginica 
#>         7.0828803         0.9045646        -1.6009717        -2.1176692 

# compare with model with no standardization of predictors
coef(lm(Sepal.Length ~ Petal.Length + Species, data = iris))
#>       (Intercept)      Petal.Length Speciesversicolor  Speciesvirginica 
#>         3.6835266         0.9045646        -1.6009717        -2.1176692 
```
