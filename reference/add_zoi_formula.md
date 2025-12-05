# Adds ZOI radii to formula

Given the structure of a formula to be fitted using a data set, the
function adds multiple radii for variables for which the zone of
influence (ZOI) or scale of effect (SOE) is to be assessed through a
statistical model. The terms representing the ZOI/SOE variables should
be identified with a pattern (see argument `pattern`). It returns the
complete formula with all variables and radii for the ZOI/SOE terms. If
`predictor_table = TRUE`, it also returns a table with information about
each predictor (e.g. whether or not it is a ZOI variable, which radius
and shape, to what type of infrastructure it corresponds), to be used by
some of the algorithms in the penalized regression modeling.

## Usage

``` r
add_zoi_formula(
  f,
  zoi_radius,
  type = "",
  pattern = "XXX",
  cumulative = "",
  separator = "_",
  remove_term = "",
  predictor_table = FALSE
)
```

## Arguments

- f:

  `[formula]`  
  A formula to be fitted in a statistical model. Here we do not need the
  actual name of each term, but all ZOI variables whose radius/scale
  shall be fitted might be added with a pattern, e.g.
  `"road_traffic_zoiXXXX"`. The pattern (e.g `"XXXX"`) needs to be set
  using the argument `pattern`.

- zoi_radius:

  `[numeric,vector]`  
  A vector of radii/scales for the zone of influence, to be added to the
  formula (e.g. `c(100, 200, 300)`).

- type:

  `[character=""]`  
  Shape(s) of the ZOI (or vector of shapes if more than one), to be
  added to the terms name (e.g. `"exp"`, `"linear"`).

- pattern:

  `[character]`  
  Pattern to be replaced in the formula, corresponding to the ZOI radii
  or ZOI shapes and radii (e.g. `"XXX"`).

- cumulative:

  `[character=""]{"cumulative", "nearest"}`  
  Default is `""`. String to be added to the ZOI terms corresponding on
  whether the variable represents the ZOI of the nearest feature
  (`cumulative = "nearest"`) or the cumulative ZOI
  (`cumulative = "cumulative"`). If `""` (default), the type of ZOI is
  taken from the variable name, already given in the formula.

- separator:

  `[character(1)="_"]`  
  Separator to be used between type and zoi_radius, when `type` is
  provided. Default is `"_"`.

- remove_term:

  `[character=""]`  
  Vector of characters with the names of the variables not to be added
  to the formula. This should be used when some specific combinations of
  variables, shapes, and radii are in principle created by the function
  but shouldbe ignored in the final formula.

- predictor_table:

  `[logical(1)=FALSE]`  
  logical. Whether or not a table should be returned with info of all
  ZOI radii, shape values, and formula terms, together with info from
  other non-ZOI predictors. This table is of special interest when
  fitting "Decay Adaptive Lasso" or other related models with
  [`fit_net_logit()`](https://ninanor.github.io/oneimpact/reference/fit_net_functions.md)
  or
  [`fit_net_clogit()`](https://ninanor.github.io/oneimpact/reference/fit_net_functions.md)
  with argument `method = "DecayAdaptiveLasso"`.

## Value

A list with both the final `formula` with all ZOI radii and shapes and a
table with the predictor information (if `predictor_table = TRUE`).

## Details

The function searches for patterns in ZOI variables as stated in the
`formula` and replaces them by combinations of strings representing the
type of ZOI (nearest, cumulative), the shape (e.g. "exp", "linear"), and
the multiple ZOI radii to be assessed. The final name of the variables
should match the names of the columns in the data set.

## Examples

``` r
# multiple radii
f <- case_ ~ strata(step_id_) + sl_*startpt_roadsXXX + sl_*startpt_cabinsXXX
add_zoi_formula(f, zoi_radius = c(1000, 2000, 3000), pattern = "XXX")
#> $formula
#> case_ ~ strata(step_id_) + sl_ * startpt_roads1000 + sl_ * startpt_roads2000 + 
#>     sl_ * startpt_roads3000 + sl_ * startpt_cabins1000 + sl_ * 
#>     startpt_cabins2000 + sl_ * startpt_cabins3000
#> <environment: 0x55976a73ec10>
#> 
#> $predictor_table
#> [1] FALSE
#> 

# multiple radii and shapes
f <- case_ ~ strata(step_id_) + sl_*startpt_roadsXXX + sl_*startpt_cabinsXXX
add_zoi_formula(f, zoi_radius = c(1000, 2000, 3000), type = "_exp_decay", pattern = "XXX")
#> $formula
#> case_ ~ strata(step_id_) + sl_ * startpt_roads_exp_decay_1000 + 
#>     sl_ * startpt_roads_exp_decay_2000 + sl_ * startpt_roads_exp_decay_3000 + 
#>     sl_ * startpt_cabins_exp_decay_1000 + sl_ * startpt_cabins_exp_decay_2000 + 
#>     sl_ * startpt_cabins_exp_decay_3000
#> <environment: 0x55976a816038>
#> 
#> $predictor_table
#> [1] FALSE
#> 
add_zoi_formula(f, zoi_radius = c(1000, 2000, 3000), type = c("_exp_decay", "_threshold"), pattern = "XXX")
#> $formula
#> case_ ~ strata(step_id_) + sl_ * startpt_roads_exp_decay_1000 + 
#>     sl_ * startpt_roads_exp_decay_2000 + sl_ * startpt_roads_exp_decay_3000 + 
#>     sl_ * startpt_roads_threshold_1000 + sl_ * startpt_roads_threshold_2000 + 
#>     sl_ * startpt_roads_threshold_3000 + sl_ * startpt_cabins_exp_decay_1000 + 
#>     sl_ * startpt_cabins_exp_decay_2000 + sl_ * startpt_cabins_exp_decay_3000 + 
#>     sl_ * startpt_cabins_threshold_1000 + sl_ * startpt_cabins_threshold_2000 + 
#>     sl_ * startpt_cabins_threshold_3000
#> <environment: 0x55976a888748>
#> 
#> $predictor_table
#> [1] FALSE
#> 

# predictor_table - adding only the nearest ZOI
f <- case_ ~ strata(step_id_) + land_use + sl_*startpt_roads_XXX + sl_*startpt_cabins_XXX
add_zoi_formula(f, zoi_radius = c(1000, 2000, 3000), type = c("exp_decay", "threshold"), pattern = "XXX",
                predictor_table = TRUE, cumulative = "nearest")
#> $formula
#> case_ ~ strata(step_id_) + land_use + sl_ * startpt_roads_exp_decay_1000 + 
#>     sl_ * startpt_roads_exp_decay_2000 + sl_ * startpt_roads_exp_decay_3000 + 
#>     sl_ * startpt_roads_threshold_1000 + sl_ * startpt_roads_threshold_2000 + 
#>     sl_ * startpt_roads_threshold_3000 + sl_ * startpt_cabins_exp_decay_1000 + 
#>     sl_ * startpt_cabins_exp_decay_2000 + sl_ * startpt_cabins_exp_decay_3000 + 
#>     sl_ * startpt_cabins_threshold_1000 + sl_ * startpt_cabins_threshold_2000 + 
#>     sl_ * startpt_cabins_threshold_3000
#> <environment: 0x55976a93c988>
#> 
#> $predictor_table
#>    is_zoi cumulative     shape zoi_radius              variable
#> 1       0       <NA>      <NA>         NA              land_use
#> 7       1    nearest exp_decay       1000  sl_ * startpt_roads_
#> 8       1    nearest exp_decay       2000  sl_ * startpt_roads_
#> 9       1    nearest exp_decay       3000  sl_ * startpt_roads_
#> 10      1    nearest threshold       1000  sl_ * startpt_roads_
#> 11      1    nearest threshold       2000  sl_ * startpt_roads_
#> 12      1    nearest threshold       3000  sl_ * startpt_roads_
#> 13      1    nearest exp_decay       1000 sl_ * startpt_cabins_
#> 14      1    nearest exp_decay       2000 sl_ * startpt_cabins_
#> 15      1    nearest exp_decay       3000 sl_ * startpt_cabins_
#> 16      1    nearest threshold       1000 sl_ * startpt_cabins_
#> 17      1    nearest threshold       2000 sl_ * startpt_cabins_
#> 18      1    nearest threshold       3000 sl_ * startpt_cabins_
#>                               term_zoi
#> 1                             land_use
#> 7   sl_ * startpt_roads_exp_decay_1000
#> 8   sl_ * startpt_roads_exp_decay_2000
#> 9   sl_ * startpt_roads_exp_decay_3000
#> 10  sl_ * startpt_roads_threshold_1000
#> 11  sl_ * startpt_roads_threshold_2000
#> 12  sl_ * startpt_roads_threshold_3000
#> 13 sl_ * startpt_cabins_exp_decay_1000
#> 14 sl_ * startpt_cabins_exp_decay_2000
#> 15 sl_ * startpt_cabins_exp_decay_3000
#> 16 sl_ * startpt_cabins_threshold_1000
#> 17 sl_ * startpt_cabins_threshold_2000
#> 18 sl_ * startpt_cabins_threshold_3000
#> 

# predictor_table 2 - adding both nearest and cumulative ZOI
f <- case_ ~ strata(step_id_) + sl_*startpt_roads_cumulative_XXX + sl_*startpt_cabins_nearest_XXX
add_zoi_formula(f, zoi_radius = c(1000, 2000, 3000), type = c("exp_decay", "threshold"), pattern = "XXX",
                predictor_table = TRUE)
#> $formula
#> case_ ~ strata(step_id_) + sl_ * startpt_roads_cumulative_exp_decay_1000 + 
#>     sl_ * startpt_roads_cumulative_exp_decay_2000 + sl_ * startpt_roads_cumulative_exp_decay_3000 + 
#>     sl_ * startpt_roads_cumulative_threshold_1000 + sl_ * startpt_roads_cumulative_threshold_2000 + 
#>     sl_ * startpt_roads_cumulative_threshold_3000 + sl_ * startpt_cabins_nearest_exp_decay_1000 + 
#>     sl_ * startpt_cabins_nearest_exp_decay_2000 + sl_ * startpt_cabins_nearest_exp_decay_3000 + 
#>     sl_ * startpt_cabins_nearest_threshold_1000 + sl_ * startpt_cabins_nearest_threshold_2000 + 
#>     sl_ * startpt_cabins_nearest_threshold_3000
#> <environment: 0x55976aaab128>
#> 
#> $predictor_table
#>    is_zoi cumulative     shape zoi_radius                        variable
#> 1       1 cumulative exp_decay       1000 sl_ * startpt_roads_cumulative_
#> 2       1 cumulative exp_decay       2000 sl_ * startpt_roads_cumulative_
#> 3       1 cumulative exp_decay       3000 sl_ * startpt_roads_cumulative_
#> 4       1 cumulative threshold       1000 sl_ * startpt_roads_cumulative_
#> 5       1 cumulative threshold       2000 sl_ * startpt_roads_cumulative_
#> 6       1 cumulative threshold       3000 sl_ * startpt_roads_cumulative_
#> 7       1    nearest exp_decay       1000   sl_ * startpt_cabins_nearest_
#> 8       1    nearest exp_decay       2000   sl_ * startpt_cabins_nearest_
#> 9       1    nearest exp_decay       3000   sl_ * startpt_cabins_nearest_
#> 10      1    nearest threshold       1000   sl_ * startpt_cabins_nearest_
#> 11      1    nearest threshold       2000   sl_ * startpt_cabins_nearest_
#> 12      1    nearest threshold       3000   sl_ * startpt_cabins_nearest_
#>                                         term_zoi
#> 1  sl_ * startpt_roads_cumulative_exp_decay_1000
#> 2  sl_ * startpt_roads_cumulative_exp_decay_2000
#> 3  sl_ * startpt_roads_cumulative_exp_decay_3000
#> 4  sl_ * startpt_roads_cumulative_threshold_1000
#> 5  sl_ * startpt_roads_cumulative_threshold_2000
#> 6  sl_ * startpt_roads_cumulative_threshold_3000
#> 7    sl_ * startpt_cabins_nearest_exp_decay_1000
#> 8    sl_ * startpt_cabins_nearest_exp_decay_2000
#> 9    sl_ * startpt_cabins_nearest_exp_decay_3000
#> 10   sl_ * startpt_cabins_nearest_threshold_1000
#> 11   sl_ * startpt_cabins_nearest_threshold_2000
#> 12   sl_ * startpt_cabins_nearest_threshold_3000
#> 

# predictor_table 3 - adding the nearest/cumulative metric within the type argument
f <- case_ ~ strata(step_id_) + sl_*startpt_roads_XXX + sl_*startpt_cabins_XXX
add_zoi_formula(f, zoi_radius = c(1000, 2000, 3000), type = c("nearest_exp_decay", "cumulative_exp_decay"), pattern = "XXX",
                predictor_table = TRUE)
#> $formula
#> case_ ~ strata(step_id_) + sl_ * startpt_roads_nearest_exp_decay_1000 + 
#>     sl_ * startpt_roads_nearest_exp_decay_2000 + sl_ * startpt_roads_nearest_exp_decay_3000 + 
#>     sl_ * startpt_roads_cumulative_exp_decay_1000 + sl_ * startpt_roads_cumulative_exp_decay_2000 + 
#>     sl_ * startpt_roads_cumulative_exp_decay_3000 + sl_ * startpt_cabins_nearest_exp_decay_1000 + 
#>     sl_ * startpt_cabins_nearest_exp_decay_2000 + sl_ * startpt_cabins_nearest_exp_decay_3000 + 
#>     sl_ * startpt_cabins_cumulative_exp_decay_1000 + sl_ * startpt_cabins_cumulative_exp_decay_2000 + 
#>     sl_ * startpt_cabins_cumulative_exp_decay_3000
#> <environment: 0x55976abf46c8>
#> 
#> $predictor_table
#>    is_zoi cumulative     shape zoi_radius              variable
#> 1       1    nearest exp_decay       1000  sl_ * startpt_roads_
#> 2       1    nearest exp_decay       2000  sl_ * startpt_roads_
#> 3       1    nearest exp_decay       3000  sl_ * startpt_roads_
#> 4       1 cumulative exp_decay       1000  sl_ * startpt_roads_
#> 5       1 cumulative exp_decay       2000  sl_ * startpt_roads_
#> 6       1 cumulative exp_decay       3000  sl_ * startpt_roads_
#> 7       1    nearest exp_decay       1000 sl_ * startpt_cabins_
#> 8       1    nearest exp_decay       2000 sl_ * startpt_cabins_
#> 9       1    nearest exp_decay       3000 sl_ * startpt_cabins_
#> 10      1 cumulative exp_decay       1000 sl_ * startpt_cabins_
#> 11      1 cumulative exp_decay       2000 sl_ * startpt_cabins_
#> 12      1 cumulative exp_decay       3000 sl_ * startpt_cabins_
#>                                          term_zoi
#> 1      sl_ * startpt_roads_nearest_exp_decay_1000
#> 2      sl_ * startpt_roads_nearest_exp_decay_2000
#> 3      sl_ * startpt_roads_nearest_exp_decay_3000
#> 4   sl_ * startpt_roads_cumulative_exp_decay_1000
#> 5   sl_ * startpt_roads_cumulative_exp_decay_2000
#> 6   sl_ * startpt_roads_cumulative_exp_decay_3000
#> 7     sl_ * startpt_cabins_nearest_exp_decay_1000
#> 8     sl_ * startpt_cabins_nearest_exp_decay_2000
#> 9     sl_ * startpt_cabins_nearest_exp_decay_3000
#> 10 sl_ * startpt_cabins_cumulative_exp_decay_1000
#> 11 sl_ * startpt_cabins_cumulative_exp_decay_2000
#> 12 sl_ * startpt_cabins_cumulative_exp_decay_3000
#> 
```
