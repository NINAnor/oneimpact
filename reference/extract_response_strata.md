# Separates elements in a statistical formula

This function separates the response variable, the strata variable, and
the covariates/predictors in a statistical formula.

## Usage

``` r
extract_response_strata(f, covars = FALSE)
```

## Arguments

- f:

  `[formula]`  
  Formula for model fitting.

- covars:

  `[logical(1)=FALSE]`  
  logical. If true, all the (other) covariates/explanatory variables in
  the formula are also appended to the output, separated from the
  response and strata variables.

## Value

A `list` with strings representing the response variable (`response`),
the stratum variable (`strata`), and possibly all explanatory variables
(`covars`), if `covars = TRUE`, as shown in the formula.

## Examples

``` r
# formula with strata
f <- formula(use ~ exp1 + exp2 + exp1:exp2 + strata(id))
extract_response_strata(f)
#> $response
#> [1] "use"
#> 
#> $strata
#> [1] "id"
#> 
extract_response_strata(f, covars = TRUE)
#> $response
#> [1] "use"
#> 
#> $strata
#> [1] "id"
#> 
#> $covars
#> [1] "exp1 + exp2 + exp1:exp2"
#> 

# formula without strata
f2 <- formula(use ~ exp1 + exp2 + exp1:exp2)
extract_response_strata(f2, covars = TRUE)
#> $response
#> [1] "use"
#> 
#> $strata
#> [1] ""
#> 
#> $covars
#> [1] "exp1 + exp2 + exp1:exp2"
#> 
```
