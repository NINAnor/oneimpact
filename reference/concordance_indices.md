# Computes concordance indices for model evaluation

These functions compute different concordance indices for
use/availability data in a train-validate-test model evaluation context.

## Usage

``` r
conditionalBoyce(
  x,
  method = c("pearson", "kendall", "spearman")[1],
  plotit = FALSE,
  errors = TRUE,
  warnings = TRUE
)

conditionalSomersD(x, errors = TRUE, warnings = TRUE)

conditionalAUC(x, errors = TRUE, warnings = TRUE)

AUC(x, errors = TRUE, warnings = TRUE)

coxnet.deviance(x, errors = TRUE, warnings = TRUE)

Cindex(x, errors = TRUE, warnings = TRUE)
```

## Arguments

- x:

  `[data.frame]`  
  A data frame with three columns: `x` (predicted values), `y`
  (use/available response, 1/0), and `strat` (stratum identifier).

- method:

  `[character(1)="pearson"]{"pearson","kendall","spearman"}`  
  Correlation method used to compute the Boyce index. Only used in
  `conditionalBoyce()`.

- plotit:

  `[logical(1)=FALSE]`  
  Whether to plot the bin frequency distribution. Only used in
  `conditionalBoyce()`.

- errors:

  `[logical(1)=TRUE]`  
  Whether to raise an error if the number of non-empty bins is
  insufficient for a reliable estimate.

- warnings:

  `[logical(1)=TRUE]`  
  Whether to emit warnings for unequal strata lengths or low bin counts.

## Value

A single numeric value representing the concordance metric:

- `conditionalBoyce()`: Pearson/Kendall/Spearman correlation of bin
  frequencies.

- `conditionalSomersD()`: Somers' D statistic (range -1 to 1).

- `conditionalAUC()`: AUC derived from Somers' D (range 0 to 1).

- `AUC()`: AUC from
  [`pROC::auc()`](https://rdrr.io/pkg/pROC/man/auc.html), ignoring
  strata.

- `coxnet.deviance()`: Cox partial deviance via
  [`glmnet::coxnet.deviance()`](https://glmnet.stanford.edu/reference/coxnet.deviance.html).

- `Cindex()`: Concordance index via
  [`glmnet::Cindex()`](https://glmnet.stanford.edu/reference/Cindex.html).

## Details

The function `conditionalAUC()` is the implementation of the AUC as
related to the Somers' D index. It accounts for strata, ideal for
conditional logistic regression, but is under testing. `AUC()` uses
[`pROC::auc()`](https://rdrr.io/pkg/pROC/man/auc.html) and does not
account for strata. `coxnet.deviance` and `Cindex` are wrappers for
[`glmnet::coxnet.deviance()`](https://glmnet.stanford.edu/reference/coxnet.deviance.html)
and
[`glmnet::Cindex()`](https://glmnet.stanford.edu/reference/Cindex.html)
using the same argument structure as the other concordance functions.
