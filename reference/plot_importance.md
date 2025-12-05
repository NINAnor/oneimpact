# Plot variable importance

Plot variable importance

## Usage

``` r
plot_importance(
  importance,
  remove_threshold = -1,
  normalize = TRUE,
  plot_type = c("barplot", "boxplot")[1],
  summ_stat = c(FALSE, "mean", "median")[1],
  by = c(FALSE, "blockH0", "variable")[1],
  scales = c("fixed", "free_y", "free_x", "free")
)
```

## Arguments

- importance:

  `[numeric]`  
  Numeric vector showing the importance score of each variable,
  resulting from
  [`variable_importance()`](https://ninanor.github.io/oneimpact/reference/variable_importance.md).

- remove_threshold:

  `[numeric(1)=0]`  
  Only importance scores above this level will be plotted. Default is 0,
  in case which all null variables will be removed. To plot all
  variables, set remove_threshold as -1 or any other negative value.

- normalize:

  `[logical(1)=TRUE]`  
  If `TRUE`, the variable importance scores are dividing the the maximum
  score, so that the score of the most important variable is set to 1.

- plot_type:

  `[character="barplot"]{"barplot", "boxplot"}`  
  Whether the plot should be a barplot (default) or a boxplot. The
  boxplot is only meaningful if `summ_stat != FALSE` (median or mean)
  and if the importance is computed for blocks (`by = "blockH0"`).

- summ_stat:

  `[string(1)=FALSE]`  
  Name of a summary statistic to be computed across blocks H0 (areas,
  herd, populations) for each variable. This is valid only when variable
  importance was computed bor each block H0, with the parameter `colH0`
  provided (i.e., not `NULL`) within the function
  [`variable_importance()`](https://ninanor.github.io/oneimpact/reference/variable_importance.md).
  Default is `FALSE`, in case which variable importance is plotted for
  each block H0. The only summary stat implemented is `"mean"` and
  `"median"` (the variable importance plotted is the average/median
  across blocks), but others might be implemented.

- by:

  `[character=FALSE]{FALSE, "blockH0", "variable"}`  
  If `by = FALSE` (default), a single plot is made. If `by = "blockH0"`,
  the variable importance if plotted for each block H0 (e.g. population,
  area) if `plot_type = "barplot"` (with the variables as a category in
  each plot), and as a single boxplot plot if `plot_type = "boxplot"`.
  If `by = "variable"`, the plot is made separately for each variable,
  with blocks H0 as a category in each plot.

- scales:

  `[character="fixed"]{"fixed", "free_x", "free_y", "free}`  
  Scale/range of the axes when facets are used in the plots, when `by`
  is not `FALSE`. See the parameter `scales` in
  [`ggplot2::facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html)
  for more information.
