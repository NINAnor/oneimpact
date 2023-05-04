#' Prediction based only on step length terms
#'
#' @export
kernel_prediction <- function(f, data,
                              kernel_vars = c("step_length", "TA"), coefs){

  all_vars <- attr(terms(f), "term.labels")
  kernel_vars <- all_vars[unlist(lapply(kernel_vars, function(x)grep(x, all_vars)))]
  kernel_vars <- kernel_vars[!grepl(":", kernel_vars)]

  f2 <- as.formula(paste0(extract_case_strata(f, other_vars=F)$case, " ~ -1+", paste0(kernel_vars, collapse = "+")))
  predVals <- model.matrix(f2, data) %*% coefs[match(kernel_vars, names(coefs))]
  return(predVals)
}
