#' Computes concordance indices for model evaluation
#'
#' These functions compute different concordance indices for use/availability data
#' in a train-validate-test model evaluation context.
#'
#' @param x `[data.frame]` \cr A data frame with three columns: `x` (predicted values),
#'   `y` (use/available response, 1/0), and `strat` (stratum identifier).
#' @param method `[character(1)="pearson"]{"pearson","kendall","spearman"}` \cr Correlation
#'   method used to compute the Boyce index. Only used in `conditionalBoyce()`.
#' @param plotit `[logical(1)=FALSE]` \cr Whether to plot the bin frequency distribution.
#'   Only used in `conditionalBoyce()`.
#' @param errors `[logical(1)=TRUE]` \cr Whether to raise an error if the number of
#'   non-empty bins is insufficient for a reliable estimate.
#' @param warnings `[logical(1)=TRUE]` \cr Whether to emit warnings for unequal strata
#'   lengths or low bin counts.
#'
#' @return A single numeric value representing the concordance metric:
#' \itemize{
#'   \item `conditionalBoyce()`: Pearson/Kendall/Spearman correlation of bin frequencies.
#'   \item `conditionalSomersD()`: Somers' D statistic (range -1 to 1).
#'   \item `conditionalAUC()`: AUC derived from Somers' D (range 0 to 1).
#'   \item `AUC()`: AUC from [pROC::auc()], ignoring strata.
#'   \item `coxnet.deviance()`: Cox partial deviance via [glmnet::coxnet.deviance()].
#'   \item `Cindex()`: Concordance index via [glmnet::Cindex()].
#' }
#'
#' @details The function [oneimpact::conditionalAUC()] is the implementation of the AUC
#' as related to the Somers' D index. It accounts for strata, ideal for conditional
#' logistic regression, but is under testing. [oneimpact::AUC()] uses [pROC::auc()] and
#' does not account for strata. `coxnet.deviance` and `Cindex` are wrappers for
#' [glmnet::coxnet.deviance()] and [glmnet::Cindex()] using the same argument structure
#' as the other concordance functions.
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
