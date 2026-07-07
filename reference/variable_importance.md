# Computes variable importance

Computation of variable importance from a bag of models. Variable
importance can be computed either by dropping variables or through
permutation of the variables values (parameter `type`), typically by
evaluating the effects on the model evaluation metric in the validation
set(s). If `type = "drop"`, each variable is dropped from the model at a
time and the variation in model evaluation metric is computed. If
`type = "permutation"`, The observations of each variable are permutated
and the variation in model evaluation metric is computed.

## Usage

``` r
variable_importance(
  x,
  data,
  samples = NULL,
  type = c("drop", "permutation")[1],
  colH0 = NULL,
  variable_block = NULL,
  n_permutations = 100,
  order = c("desc", "asc", FALSE)[1],
  metric = NULL,
  plot = FALSE,
  ss = 1,
  remove_threshold = 0
)
```

## Arguments

- x:

  `[list]`  
  Bag of models, result of
  [`bag_models()`](https://ninanor.github.io/oneimpact/reference/bag_models.md).
  It contains multiple information from the models, such as formula,
  weights, coefficients, and metric used to evaluate the models.

- data:

  `[data.frame]`  
  Complete data set to which the models were applied.

- samples:

  `[list]`  
  List of samples used to fit the models in the bag. The list contains
  at least three elements: train, test, and validate. Each elements
  might have several elements, each representing the lines of `data` to
  be sampled for each resample. Typically, this is computed by the
  function
  [`create_resamples()`](https://ninanor.github.io/oneimpact/reference/create_resamples.md).

- type:

  `[character(1)="drop"]{"drop", "permutation"}`  
  Type of computation for variable importance. If `type = "drop"`
  (default), each variable is dropped from the model at a time and the
  variation in model evaluation metric is computed. If
  `type = "permutation"`, the observations of each variable are
  permutated and the variation in the model evaluation metric is
  computed.

- colH0:

  `[string(1)=NULL]`  
  String with the name of the column in `data` representing the blockH0,
  in case we want the variable importance to be evaluated for each
  block. Default is `NULL`, in case variable importance is assessed for
  all the data.

- variable_block:

  `[character,NULL]`  
  Optional grouping variable for computing importance at the block level
  rather than per-variable. Default is `NULL`.

- n_permutations:

  `[numeric(1)=100]`  
  Number of permutations, if `type = "permutation"`.

- order:

  `[character,logical(1)="desc"]{"desc", "asc", FALSE}`  
  Whether or not to order the output variables according to descending
  (`order = "desc"`) or ascending order of variable importance
  (`order = "asc"`). If `FALSE`, the variables are shown in the same
  order as present in the bag of models, `x`.

- metric:

  `[function,character=NULL]`  
  Metric function used to evaluate model performance. If `NULL`
  (default), the metric stored in the bag `x` is used.

- plot:

  `[logical(1)=FALSE]`  
  Should variable importance be plotted? Default is `FALSE`.

- ss:

  `[integer(1)=1]`  
  Index of the sample (resample) to use for validation when `samples` is
  provided. Default is `1`.

- remove_threshold:

  `[numeric(1)]`  
  Threshold for excluding variable with little importance in the
  variable importance plot (i.e. only considered if `plot = TRUE`). See
  more in
  [`plot_importance()`](https://ninanor.github.io/oneimpact/reference/plot_importance.md).

## Value

A data frame with variable names and their importance scores, optionally
ordered by importance and/or plotted.

## See also

For plotting variable importance, see
[`plot_importance()`](https://ninanor.github.io/oneimpact/reference/plot_importance.md).
