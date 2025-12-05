# Explore potential hierarchical blocks before sampling or spatial stratification

Function to explore the number of cases and observations for the
different sampling units possibly used as the base H0 hierarchical
level, such as population ID, study area, animal ID, or year, before
spatial stratification or creating samples for the bootstrapped
approach. The function can help understand how imbalanced is the data
across H0 levels used for validation.

## Usage

``` r
explore_blocks_pre(data, colH0, animal_id = NULL, col_case = "case")
```

## Arguments

- data:

  `[data.frame,tibble]`  
  Complete data set to be analyzed.

- colH0:

  `[character]`  
  Name of the column in `data` to be used as the H0 hierarchical level,
  intended for model validation.

- animal_id:

  `[character]`  
  Name of the column in `data` representing animal ID. If `NULL`
  (default), summaries are not created for individuals.

- col_case:

  `[string(1)="case"]`  
  Name of the column in `data` representing the case or used/available
  points. Default is `"case"`.

## Examples

``` r
# read data
data("reindeer_ssf")

# explore blocks - animal ID as block H0
explore_blocks_pre(reindeer_ssf, "original_animal_id", col_case = "case_")
#> # A tibble: 9 × 3
#>   original_animal_id     n n_presences
#>                <dbl> <int>       <int>
#> 1               3358  5346         486
#> 2               3361  5412         492
#> 3               3362  2706         246
#> 4               3364  5214         474
#> 5               3372  2706         246
#> 6               3378  2706         246
#> 7               6331  2706         246
#> 8               6333  2233         203
#> 9               6335  2706         246

# explore blocks - year as block H0
library(lubridate)
#> 
#> Attaching package: ‘lubridate’
#> The following objects are masked from ‘package:terra’:
#> 
#>     intersect, union
#> The following objects are masked from ‘package:base’:
#> 
#>     date, intersect, setdiff, union
reindeer_ssf |>
  dplyr::mutate(year = lubridate::year(t1_)) |>
  explore_blocks_pre("year", col_case = "case_")
#> # A tibble: 3 × 3
#>    year     n n_presences
#>   <dbl> <int>       <int>
#> 1  2007 16071        1461
#> 2  2008  8019         729
#> 3  2009  7645         695

# year as block H0 + animal ID
reindeer_ssf |>
  dplyr::mutate(year = lubridate::year(t1_)) |>
  explore_blocks_pre("year", animal_id = "original_animal_id", col_case = "case_")
#> # A tibble: 12 × 4
#>     year original_animal_id     n n_presences
#>    <dbl>              <dbl> <int>       <int>
#>  1  2007               3358  2640         240
#>  2  2007               3361  2706         246
#>  3  2007               3362  2706         246
#>  4  2007               3364  2607         237
#>  5  2007               3372  2706         246
#>  6  2007               3378  2706         246
#>  7  2008               3358  2706         246
#>  8  2008               3361  2706         246
#>  9  2008               3364  2607         237
#> 10  2009               6331  2706         246
#> 11  2009               6333  2233         203
#> 12  2009               6335  2706         246
```
