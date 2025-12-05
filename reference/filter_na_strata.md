# Remove missing values and ensure case-control per stratum

Remove strata for which there are only 1s (used points, presence) or
only 0s (available points, background or absence).

## Usage

``` r
filter_na_strata(f, data)
```

## Arguments

- f:

  `[formula]`  
  Formula for model fitting.

- data:

  `[data.frame]`  
  Data set to be checked.

## Value

Cleaned `data.frame`, removing from the input `data` the rows with NA in
any of the columns and all the strata for which there are only presences
or only absences.

## References

[`amt::remove_incomplete_strata()`](https://rdrr.io/pkg/amt/man/remove_incomplete_strata.html)

## Examples

``` r
library(survival)
# create test dataset
test1 <- data.frame(case = c(0,0,1,0,1,1,0,1,NA,1),
                    x = c(0,2,1,2,1,0,1,0,0,1),
                    step_id = c(0,0,0,0,1,1,1,2,2,2))
# Remove NAs
f <- case ~ x + strata(step_id)
filter_na_strata(f, test1)
#>   case x step_id
#> 1    0 0       0
#> 2    0 2       0
#> 3    1 1       0
#> 4    0 2       0
#> 5    1 1       1
#> 6    1 0       1
#> 7    0 1       1

# Fit a stratified model
if(FALSE) {
  # with no NAs; necessary if using glmnet
  coxph(Surv(rep(1, length(case)), case) ~ x + strata(step_id), filter_na_strata(f, test1))
  # This differs from that
  coxph(Surv(rep(1, length(case)), case) ~ x + strata(step_id), test1)
}
```
