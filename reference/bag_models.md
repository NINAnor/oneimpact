# Summary of a bag of models

Summary of a bag of models

## Usage

``` r
bag_models(
  fitted,
  data,
  score2weight = score2weight_min_invmean,
  metric = fitted$metric,
  weights_col = c("validation_score", "habitat_validation_score")[1],
  score_threshold = 0.7,
  weights_function = NULL,
  out_dir
)

w_strech_maxmin_squared(x)

w_strech_max_squared(x)

score2weight_mean(x, col = "validation_score", score_threshold = 0.7)

score2weight_min_mean(x, col = "validation_score", score_threshold = 0.7)

score2weight_invmean(x, col = "validation_score", score_threshold = 0.7)

score2weight_min_invmean(x, col = "validation_score", score_threshold = 0.7)
```

## Arguments

- data:

  `[data.frame,tibble]`  
  Complete data set used to fit the models.

- score2weight:

  Function to set validation scores into weights, with two arguments: x,
  the result of one model of the bag, and col, the column to be used for
  setting the scores. See the argument `weights_col`.

- weights_col:

  Column to use for scores. One of `"validation_score"` or
  `"habitat_validation_score"`.
