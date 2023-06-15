#' Prediction based only on step length terms
#'
#' @export
kernel_prediction <- function(f, data,
                              kernel_vars = c("step_length", "ta"),
                              coefs){

  # get movement/kernel variables from the formula
  # ignoring the interactions
  all_vars <- attr(terms(f), "term.labels")
  kernel_variables <- all_vars[grepl(paste(kernel_vars, collapse = "|"), all_vars)]
  kernel_variables <- kernel_variables[!grepl(":", kernel_variables)]

  # returning prediction of only this variables, based on the fitted coefficients
  f2 <- as.formula(paste0(extract_response_strata(f, other_vars = F)$response, " ~ -1 + ",
                          paste0(kernel_variables, collapse = "+")))
  pred_vals_kernel <- model.matrix(f2, data) %*% coefs[match(kernel_variables, names(coefs))]
  return(pred_vals_kernel)
}
