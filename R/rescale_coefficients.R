#' Rescale standardized coefficients back to their original range after model fitting
#'
#' Predictor variables are often standardized to be included in statistical models
#' and allow comparison of the effect sizes for different predictors. This functions
#' scales the fitted models coefficients back to the original scale of the predictors,
#' to allow ecological interpretation.
#'
#' @param model \cr Fitted model, i.e. the object created by a model fit funtion
#' such as "lm", "glm", or "coxph".
#' @param data \cr The original `data.frame` with the data used to fit the model.
#' @param bag \cr A bag of models, as result of the [oneimpact::bag_models()] function.
#' @param standardize `[logical(1)=TRUE]` If `TRUE` (Default), the coefficients are standardized.
#' If `FALSE`, the coefficients are standardized. Only numeric coefficients are
#' standardized.
#'
#' @returns A vector of rescaled coefficients for the input model.
#'
#' @name rescale_coefficients
#'
#' @examples
#' library(dplyr)
#'
#' # standardize predictors
#' iris_std <- iris |>
#'   dplyr::mutate(across(2:4, ~ scale(.x)))
#' # fit model
#' m1 <- lm(Sepal.Length ~ Petal.Length + Species, data = iris_std)
#' summary(m1)
#'
#' # rescale coefficients
#' (resc_cf <- rescale_coefficients(m1, iris))
#'
#' # compare with model with no standardization of predictors
#' coef(lm(Sepal.Length ~ Petal.Length + Species, data = iris))
#'
#' @export
# rescale_coefficients <- function(model, data) {
#
#   if("glm" %in% class(model)) {
#     beta <- rescale_coefs.glm(model, data)
#   } else {
#     if("lm" %in% class(model))
#       beta <- rescale_coefs.lm(model, data)
#   }
#
#   if("coxph" %in% class(model))
#     beta <- rescale_coefs.coxph(model, data)
#
#   beta
# }
rescale_coefficients <- function(...) {
  UseMethod("rescale_coefficients")
}

#' @rdname rescale_coefficients
#' @export
rescale_coefficients.coxph <- function(model, data, ...) {

  # model inputs
  # m_vars <- all.vars(model$formula)[-1] # remove reponse variable
  # model matrix with data
  M <- stats::model.matrix(model, data)

  beta <- stats::coef(model) ## inherit names etc.

  # get sigma from the model matrix
  factors <- names(attributes(M)$contrasts)
  vars <- attributes(M)$dimnames[[2]]
  sigma <- apply(M, 2, stats::sd)
  sigma <- ifelse(grepl(factors, vars), 1, sigma)

  # rescale
  beta <- beta/sigma

  # beta2[1]  <- sigma[1]*beta[1]+mu[1]-sum(beta2[-1]*mu[-1])
  beta
}

#' @rdname rescale_coefficients
#' @export
rescale_coefficients.lm <- function(model, data, ...) {

  # model inputs
  m_vars <- all.vars(attr(model$terms, "variables"))[-1] # remove reponse variable
  # data - only the interest variables
  dd <- data[m_vars]

  # model matrix with data
  M <- stats::model.matrix(model, data)

  beta <- stats::coef(model) ## inherit names etc.

  # get sigma from the model matrix
  factors <- names(attributes(M)$contrasts)
  terms <- attributes(M)$dimnames[[2]]
  sigma <- apply(M, 2, sd, na.rm = T)
  sigma_rawdata <- apply(subset(dd, select = -get(factors)), 2, stats::sd, na.rm = T)
  sigma <- ifelse(grepl(paste("(Intercept)", factors, sep = "|"), terms), 1, NA)
  sigma[is.na(sigma)] <- sigma_rawdata

  # rescale
  beta <- beta/sigma

  # beta2[1]  <- sigma[1]*beta[1]+mu[1]-sum(beta2[-1]*mu[-1])
  beta
}

#' @rdname rescale_coefficients
#' @export
rescale_coefficients.glm <- function(model, data, ...) {

  # model inputs
  m_vars <- all.vars(model$formula)[-1] # remove reponse variable
  # data - only the interest variables
  dd <- data[m_vars]

  # model matrix with data
  M <- stats::model.matrix(model, data)

  beta <- stats::coef(model) ## inherit names etc.

  # get sigma from the model matrix
  factors <- names(attributes(M)$contrasts)
  terms <- attributes(M)$dimnames[[2]]
  sigma <- apply(M, 2, sd, na.rm = T)
  sigma_rawdata <- apply(subset(dd, select = -get(factors)), 2, stats::sd, na.rm = T)
  sigma <- ifelse(grepl(paste("(Intercept)", factors, sep = "|"), terms), 1, NA)
  sigma[is.na(sigma)] <- sigma_rawdata

  # rescale
  beta <- beta/sigma

  # beta2[1]  <- sigma[1]*beta[1]+mu[1]-sum(beta2[-1]*mu[-1])
  beta
}

#' @rdname rescale_coefficients
#' @export
rescale_coefficients.bag  <- function(bag, data, tostd = TRUE, ...) {

  # covariates
  m_covars <- all.vars(bag$formula_no_strata, unique = F)[-1]
  # are they numeric?
  repeated <- m_covars[which(duplicated(m_covars))]
  rep_times <- ifelse(names(bag$numeric_covs) %in% repeated, 2, 1) ## CORRECT IF THERE ARE MORE THAN TWO TERMS WITH THE SAME VARIABLE
  numeric_covs <- rep(bag$numeric_covs, times = rep_times)

  # model matrix with data
  M <- stats::model.matrix(bag$formula_no_strata, data)

  # variables and terms
  terms_order <- attributes(M)$assign
  terms_order <- terms_order[terms_order > 0]
  vars_formula <- rep(m_covars, times = unname(table(terms_order)))
  numeric_vars_order <- rep(numeric_covs, times = unname(table(terms_order)))
  # numeric variables in the formula
  m_covars_num <- m_covars[numeric_covs]
  vars_formula_num <- vars_formula[vars_formula %in% m_covars_num]

  # terms_order <- table(m_covars)
  # terms_order <- terms_order[order(match(names(terms_order), m_covars))]

  # coefficients
  coef <- bag$coef

  # SDs
  sds <- bag$data_summary
  sds <- sds[rownames(sds) == "sd", colnames(sds) %in% m_covars]
  # sds_all <- unlist(rep(sds, unname(table(terms_order)))) |>
  #   as.numeric()
  sds_all <- sds[match(vars_formula, colnames(sds))]
  # sds_all <- unlist(rep(sds, terms_order)) |>
  #   as.numeric()
  sds_all[numeric_vars_order == FALSE] <- 1
  names(sds_all) <- vars_formula
  sds_all <- unlist(sds_all)

  # standardized coefs
  if(tostd)
    new_coef <- to_std(coef, sds_all) else
      new_coef <- from_std(coef, sds_all)

  new_coef
}

# function to standardize coefficients that are not standardized
# here coef and sd should be the same length
#' @export
to_std <- function(coef, sd) {
  coef * sd
}

#' @export
from_std <- function(coef, sd) {
  coef / sd
}
