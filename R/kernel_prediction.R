#' Prediction based only on step length terms
#'
#' Computes model predictions based exclusively on movement kernel variables
#' (e.g. step length and turning angle), using the fitted coefficients from
#' a conditional logistic regression model.
#'
#' @param f `[formula]` \cr Model formula, used to identify kernel variable terms.
#' @param data `[data.frame]` \cr Data set on which predictions are made.
#' @param kernel_vars `[character=c("step_length","ta")]` \cr Names of the movement
#'   kernel variables in the model.
#' @param coefs `[numeric]` \cr Named numeric vector of fitted model coefficients.
#'
#' @return A numeric matrix of predicted values based on kernel variables only.
#'
#' @seealso [oneimpact::fit_net_clogit()], [oneimpact::bag_fit_net_clogit()]
#'
#' @export
kernel_prediction <- function(f, data,
                              kernel_vars = c("step_length", "ta"),
                              coefs) {

  # get movement/kernel variables from the formula
  # ignoring the interactions
  all_vars <- attr(terms(f), "term.labels")
  kernel_variables <- all_vars[grepl(paste(kernel_vars, collapse = "|"), all_vars)]
  kernel_variables <- kernel_variables[!grepl(":", kernel_variables)]
  # make sure to remove strata
  kernel_variables <- kernel_variables[grep("strata", kernel_variables, invert = TRUE)]

  # returning prediction of only this variables, based on the fitted coefficients
  f2 <- as.formula(paste0(extract_response_strata(f, covars = F)$response, " ~ -1 + ",
                          paste0(kernel_variables, collapse = "+")))
  pred_vals_kernel <- model.matrix(f2, data) %*% coefs[match(kernel_variables, names(coefs))]
  return(pred_vals_kernel)
}
