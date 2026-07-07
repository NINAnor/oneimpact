# Plot the coefficients of bags of models

Visualizes the distribution of coefficients across all models in a bag,
optionally weighted, filtered by term, and plotted as bars, points, or
histograms. Supports standardized, unstandardized, or raw coefficients.

## Usage

``` r
plot_coef(
  bag,
  terms = "all",
  models = 1:bag$n_above_threshold,
  weighted = TRUE,
  remove_weight_zero = TRUE,
  what = c("all_models", "average")[1],
  plot_type = c("bars", "points", "histogram")[1],
  remove_low = -1,
  remove_high = Inf,
  std = c(FALSE, "std", "unstd")[1],
  data = NULL,
  order_zoi_radius = FALSE,
  show_legend = FALSE
)
```

## Arguments

- bag:

  `[bag,list]`  
  Bag of models as returned by
  [`bag_models()`](https://ninanor.github.io/oneimpact/reference/bag_models.md).

- terms:

  `[character="all"]`  
  Which terms to include in the plot. Can be `"all"` or a pattern.

- models:

  `[integer]`  
  Model indices to include. Default is `1:bag$n_above_threshold`.

- weighted:

  `[logical(1)=TRUE]`  
  If `TRUE`, coefficients are weighted by model weights.

- remove_weight_zero:

  `[logical(1)=TRUE]`  
  If `TRUE`, models with zero weight are excluded when weighting.

- what:

  `[character="all_models"]{"all_models", "average"}`  
  Plot coefficients for all models or the weighted average.

- plot_type:

  `[character="bars"]{"bars", "points", "histogram"}`  
  Plot type to use.

- remove_low:

  `[numeric(1)=-1]`  
  Minimum coefficient value to include.

- remove_high:

  `[numeric(1)=Inf]`  
  Maximum coefficient value to include.

- std:

  `[logical,character="FALSE"]{FALSE,"std","unstd"}`  
  If `"std"`, use standardized coefficients; if `"unstd"`,
  unstandardized; if `FALSE` (default), keep raw coefficients.

- data:

  `[data.frame]`  
  Data used to standardize or unstandardize coefficients when `std` is
  not `FALSE`.

- order_zoi_radius:

  `[logical(1)=FALSE]`  
  If `TRUE`, order ZOI terms by radius before plotting.

- show_legend:

  `[logical(1)=FALSE]`  
  Whether to show a legend in the plot.

## Value

A ggplot2 object representing the coefficient plot.
