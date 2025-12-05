# Computes the conditional Boyce index for model evaluation

These functions compute different concordance indices.

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

  A data.frame with three columns: x, the predicted values; y, the case
  variable (use vs. available, 1/0); and strat, the stratum.

## Details

The function `conditionalAUC()` is the implementation of the computation
of the Area Under the Curve as related to the Sommers'D index, as
described in
https://cran.r-project.org/web/packages/survival/vignettes/concordance.pdf.
It is implemented by accounting for strata, ideal for conditional
logistic regression, but it is under testing. The function `AUC()` uses
the [`pROC::auc()`](https://rdrr.io/pkg/pROC/man/auc.html) function and
does not account for strata. Functions conditionalBoyce and
conditionalSomersD both account for strata, but are under checking. The
functions coxnet.deviance and Cindex are wrappers for the same functions
from glmnet, which are suitable for cox models, but here they use the
same arguments as all the other concordance variables.

## References

https://cran.r-project.org/web/packages/survival/vignettes/concordance.pdf

## See also

conditionalBoyce
