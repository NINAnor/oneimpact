#' Fits an issf using glmnet
#'
#' @param f formula
#' @param data data
#'
#' @export
net_issf <- function(f, data, ...) {

  # extract strata variable and other columns
  wcols <- extract_case_strata(f, other_vars = TRUE)
  # response values
  p <- data[[wcols$case]]
  # strata values
  strata <- data[[wcols$strata]]

  if (anyNA(data)) stop("NA values in data table. Please remove them and rerun.")
  lp <- split(p, strata)
  test <- range(unlist(lapply(lp, mean)))
  if (anyNA(test)) stop("NA values in the response. Please remove them and rerun.")
  if (test[1]==0) stop("There are strata without event. Please remove them and rerun.")
  if (test[1]==1) stop("There are strata without control. Please remove them and rerun.")

  fit <- glmnet::glmnet(model.matrix(as.formula(paste0(wcols$case, " ~ -1+", wcols$other_vars)), data),
                        glmnet::stratifySurv(survival::Surv(rep(1, length(p)), p), strata),
                        family = "cox", ...)
  return(fit)
}
