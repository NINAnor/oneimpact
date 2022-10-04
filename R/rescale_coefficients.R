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
#'
#' @returns A vector of rescaled coefficients for the input model.
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
#' resc_cf <- rescale_coefficients(m1, iris)
#'
#' # compare with model with no standardization of predictors
#' coef(lm(Sepal.Length ~ Petal.Length + Species, data = iris))
#'
#' @export
rescale_coefficients <- function(model, data) {

  if("glm" %in% class(model))
    beta <- rescale_coefs.glm(model, data)

  if("lm" %in% class(model))
    beta <- rescale_coefs.lm(model, data)

  if("coxph" %in% class(model))
    beta <- rescale_coefs.coxph(model, data)

  beta
}
# rescale_coefficients <- function(model, data) {
#   UseMethod("rescale_coefs")
# }

rescale_coefs.coxph <- function(model, data) {

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

rescale_coefs.lm <- function(model, data) {

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

rescale_coefs.glm <- function(model, data) {

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


# model = multifits$fits[[1]]$models[[6]]
# data = dat
# beta = coef(model)
