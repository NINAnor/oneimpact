#' Rescale standardized coefficients back to their original range after model fitting
#'
#' @export
rescale_coefficients <- function(model, data) {

  if("glm" %in% class(model))
    beta <- rescale_coefs.glm(model, data)

  # if("lm" %in% class(model))
  #   beta <- rescale_coefs.lm(model, data)

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
  M <- model.matrix(model, data)

  beta <- coef(model) ## inherit names etc.

  # get sigma from the model matrix
  factors <- names(attributes(M)$contrasts)
  vars <- attributes(M)$dimnames[[2]]
  sigma <- apply(M, 2, sd)
  sigma <- ifelse(grepl(factors, vars), 1, sigma)

  # rescale
  beta <- beta/sigma

  # beta2[1]  <- sigma[1]*beta[1]+mu[1]-sum(beta2[-1]*mu[-1])
  beta
}

rescale_coefs.lm <- function(model, data) {

  # model inputs
  m_vars <- all.vars(model$formula)[-1] # remove reponse variable
  # data - only the interest variables
  dd <- data[m_vars]

  # model matrix with data
  M <- model.matrix(model, data)

  beta <- coef(model) ## inherit names etc.

  # get sigma from the model matrix
  factors <- names(attributes(M)$contrasts)
  terms <- attributes(M)$dimnames[[2]]
  sigma <- apply(M, 2, sd, na.rm = T)
  sigma_rawdata <- apply(subset(dd, select = -get(factors)), 2, sd, na.rm = T)
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
  M <- model.matrix(model, data)

  beta <- coef(model) ## inherit names etc.

  # get sigma from the model matrix
  factors <- names(attributes(M)$contrasts)
  terms <- attributes(M)$dimnames[[2]]
  sigma <- apply(M, 2, sd, na.rm = T)
  sigma_rawdata <- apply(subset(dd, select = -get(factors)), 2, sd, na.rm = T)
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
