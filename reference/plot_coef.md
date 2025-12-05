# Plot the coefficients of bags of models

Options - one term, multiple bars side by side Colors for negative and
positive

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
