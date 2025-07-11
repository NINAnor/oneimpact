% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_net_clogit.R, R/fit_net_logit.R
\name{bag_fit_net_functions}
\alias{bag_fit_net_functions}
\alias{bag_fit_net_clogit}
\alias{bag_fit_net_logit}
\title{Fit a bag of conditional logistic regression/SSF/iSSF models with penalized regression in a train-validate-test setup}
\usage{
bag_fit_net_clogit(
  f,
  data,
  samples,
  subset_samples = 1:length(samples$train),
  kernel_vars = c("step_length", "ta"),
  metric = c("coxnet.deviance", "Cindex", "conditionalAUC", "conditionalSomersD")[1],
  metrics_evaluate = c("coxnet.deviance", "Cindex", "conditionalAUC"),
  method = c("Lasso", "Ridge", "AdaptiveLasso", "DistanceDecayLasso", "DDLasso",
    "OneZOI-AdaptiveLasso", "OZ-AdaptiveLasso", "Grouped-AdaptiveLasso",
    "G-AdaptiveLasso", "HypothesisDriven-AdaptiveLasso", "HD-AdaptiveLasso",
    "ElasticNet")[1],
  standardize = c("internal", "external", FALSE)[1],
  alpha = NULL,
  penalty.factor = NULL,
  predictor_table = NULL,
  na.action = "na.pass",
  out_dir_file = NULL,
  parallel = c(FALSE, "foreach", "mclapply")[1],
  mc.cores = 2L,
  verbose = FALSE,
  ...
)

bag_fit_net_logit(
  f,
  data,
  samples,
  metric = c(AUC, conditionalBoyce, conditionalSomersD, conditionalAUC)[[1]],
  method = c("Lasso", "Ridge", "AdaptiveLasso", "DistanceDecayLasso", "ElasticNet")[1],
  standardize = c("internal", "external", FALSE)[1],
  alpha = NULL,
  penalty.factor = NULL,
  predictor_table = NULL,
  na.action = "na.pass",
  out_dir_file = NULL,
  parallel = c(FALSE, "foreach", "mclapply")[1],
  mc.cores = 2L,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{standardize}{internal = internal glmnet standaridization, i.e. using glmnet with argument standardize = TRUE.
This also standardizes dummy variables, but returns the estimated coefficients back to the original scale.
This however can cause baises in the estimates because of the bias-variance tradeoff that L1 and L1 regularization
methods try to minimize.
See more info in https://stackoverflow.com/questions/17887747/how-does-glmnets-standardize-argument-handle-dummy-variables
external = glmnet is called with argument standardize = FALSE, but standization is done by the
bag_fit_net_logit function. Return coefs in the original scale?? Implement.
If FALSE, no standardization of predictors is done.}

\item{mc.cores}{Only relevant if \code{parallel == "mclapply"}. If \code{parallel == "foreach"}, cores must
be assigned before running \code{fit_multi_net_logit()} using \code{\link[parallel:makeCluster]{parallel::makeCluster()}} and
\code{\link[doParallel:registerDoParallel]{doParallel::registerDoParallel()}}.}

\item{...}{Options for net_logit and glmnet}

\item{subset}{\verb{[vector]} \cr Vector of samples to be run (e.g. \code{c(1,2,4,5)} or \code{3:10}).
By default, all the samples in \code{samples} are run.}
}
\description{
Fit a bag of conditional logistic regression/SSF/iSSF models with penalized regression in a train-validate-test setup

Fit a bag of logistic regression/RSF models with penalized regression in a train-validate-test setup
}
