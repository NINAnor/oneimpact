# Helper functions for weight calculation in bag models

These functions compute weights for ensemble model combinations based on
validation scores. They handle score thresholding and weight
normalization differently:

- `score2weight_*` functions convert validation scores to weights

- `w_strech_*` functions normalize weights after they are computed

## Usage

``` r
w_strech_maxmin_squared(x)

w_strech_max_squared(x)

score2weight_mean(x, col = "validation_score", score_threshold = 0.7)

score2weight_min_mean(x, col = "validation_score", score_threshold = 0.7)

score2weight_invmean(x, col = "validation_score", score_threshold = 0.7)

score2weight_min_invmean(x, col = "validation_score", score_threshold = 0.7)
```

## Arguments

- x:

  `[list or numeric]`  
  For `score2weight_*` functions: a model result list with a column
  specified by `col` parameter. For `w_strech_*` functions: a numeric
  vector of weights.

- col:

  `[character(1)="validation_score"]`  
  Column name in `x` containing validation scores.

- score_threshold:

  `[numeric(1)=0.7]`  
  Minimum validation score threshold. Models with any validation score
  below this are assigned zero weight.

## Value

Numeric vector or value representing normalized weights or summary
statistics.

## Details

**Score to weight conversion functions:**

- `score2weight_mean`: Uses mean of all validation scores

- `score2weight_min_mean`: Uses mean only if all scores above threshold

- `score2weight_invmean`: Uses inverse of mean (harmonic mean-related)

- `score2weight_min_invmean`: Uses inverse mean only if all scores above
  threshold

**Weight stretching/normalization functions:**

- `w_strech_maxmin_squared`: Rescales to 0,1, squares, then normalizes
  to sum to 1

- `w_strech_max_squared`: Divides by max, squares, then normalizes to
  sum to 1

## See also

[`bag_models()`](https://ninanor.github.io/oneimpact/reference/bag_models.md),
[`data_summary()`](https://ninanor.github.io/oneimpact/reference/data_summary.md)
