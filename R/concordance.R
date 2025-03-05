#' Computes the conditional Boyce index for model evaluation
#'
#' @param x A data.frame with three columns: x, the predicted values; y, the case variable
#' (use vs. available, 1/0); and strat, the stratum.
#'
#' @details The function [oneimpact::conditionalAUC()] is the implementation of the computation of the Area Under the Curve as related
#' to the Sommers'D index, as described in https://cran.r-project.org/web/packages/survival/vignettes/concordance.pdf.
#' It is implemented by accounting for strata, ideal for conditional logistic regression, but it is under testing.
#' The function [oneimpact::AUC()] uses the [pROC::auc()] function and does not account for strata.
#' Functions conditionalBoyce and conditionalSomersD both account for strata, but are under checking.
#' The functions coxnet.deviance and Cindex are wrappers for the same functions from glmnet, which are
#' suitable for cox models, but here they use the same arguments as all the other concordance variables.
#'
#' @name concordance_indices
#' @export
conditionalBoyce <- function(x,
                             method=c("pearson", "kendall", "spearman")[1],
                             plotit = FALSE,
                             errors = TRUE,
                             warnings = TRUE){
  # x: dataframe with x: predicted value, y:use vs. available, strat: stratum
  # https://cran.r-project.org/web/packages/survival/vignettes/concordance.pdf

  # split by strata
  lx <- split(x, x$strat)
  # length of each strata
  rng <- range(unlist(lapply(lx, nrow)))

  # test for problems in sampling
  test <- as.data.frame(table(unlist(lapply(lx, function(x){ which(x$y[order(x$x)]==1)}))))
  if (errors) {
    if (nrow(test) < 3) stop(paste0("Only ", nrow(test), " bins contain observations, the correlation is unreliable."))
  }
  if (warnings){
    if (rng[1] != rng[2]) warning("Not all strata have the same length.")
    if (length(lx) > 1 & nrow(test) < (rng[1]-1)) warning(paste0("Only ", nrow(test), " bins contain observations, the correlation is likely not reliable."))
  }

  # compute correlation
  test2 <- data.frame(bins = c(1:rng[2]), frequency = 0)
  test2$frequency[match(as.character(test[,1]), as.character(test2$bins))] <- test$Freq
  b <- cor(test2$bins, test2$frequency, method = method)
  if (plotit) plot(test2)
  return(b)
}

#' Computes concordance indices for model evaluation
#'
#' These functions compute different concordance indices.
#'
#' @references https://cran.r-project.org/web/packages/survival/vignettes/concordance.pdf
#' @seealso [oneimpact::conditionalBoyce]
#'
#' @rdname concordance_indices
#' @export
conditionalSomersD <- function(x,
                               errors = TRUE,
                               warnings = TRUE){
  #x: dataframe with x: predicted value, y:use vs. available, strat: stratum
  #https://cran.r-project.org/web/packages/survival/vignettes/concordance.pdf
  lx <- split(x, x$strat)
  C  <- sum(unlist(lapply(lx, function(x){sum(x$x[x$y==1]>x$x[x$y==0])})))
  D  <- sum(unlist(lapply(lx, function(x){sum(x$x[x$y==1]<x$x[x$y==0])})))
  Tx  <- sum(unlist(lapply(lx, function(x){sum(x$x[x$y==1]==x$x[x$y==0])})))
  d <- (C-D)/(C+D+Tx)
  return(d)
}

#' @rdname concordance_indices
#' @export
conditionalAUC <- function(x,
                           errors = TRUE,
                           warnings = TRUE) {
  d <- conditionalSomersD(x)
  auc <- (d + 1)/2
  auc
}

# Another possible function for presence/absence responses
# Default cost for binomial outcomes in boot::cv.glm
cost <- function(r, pi = 0, na.rm = TRUE) mean(abs(r-pi) > 0.5, na.rm = na.rm)

# Different version of the AUC, using an external package
#' @rdname concordance_indices
#' @export
AUC <- function(x,
                errors = TRUE,
                warnings = TRUE) {
  auc_val <- as.numeric(pROC::auc(pROC::roc(x$y, x$x, quiet = TRUE)))
  auc_val
}

# the same as glmnet::coxnet.deviance, but using the same input/parameters as the other functions
#' @rdname concordance_indices
#' @export
coxnet.deviance <- function(x,
                            errors = TRUE,
                            warnings = TRUE) {
  yy <- glmnet::stratifySurv(survival::Surv(rep(1, length(x$y)), x$y), x$strat)
  dev_val <- glmnet::coxnet.deviance(pred = x$x, y = yy)
  dev_val
}

# the same as glmnet::Cindex, but using the same input/parameters as the other functions
#' @rdname concordance_indices
#' @export
Cindex <- function(x,
                   errors = TRUE,
                   warnings = TRUE) {
  yy <- glmnet::stratifySurv(survival::Surv(rep(1, length(x$y)), x$y), x$strat)
  dev_val <- glmnet::Cindex(pred = x$x, y = yy)
  dev_val
}
