# Compute data summary statistics

Utility functions to compute summary statistics for numeric and
character variables. Used internally by
[`bag_models()`](https://ninanor.github.io/oneimpact/reference/bag_models.md)
to create data summaries.

## Usage

``` r
data_summary(x, na.rm = TRUE)

data_summary_char(x)
```

## Arguments

- x:

  `[numeric or character vector]`  
  Vector of values to summarize.

- na.rm:

  `[logical(1)=TRUE]`  
  Should NA values be removed before computation?

## Value

Named numeric vector with summary statistics.

## Details

`data_summary()` computes min, quantiles (1%, 2.5%, 25%, 50%, 75%,
97.5%, 99%), max, mean, and standard deviation for numeric variables.

`data_summary_char()` computes mode and quantile-like positions for
character/factor variables, returning results with compatible naming.
