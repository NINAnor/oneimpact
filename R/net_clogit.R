#' Fits a conditional logistic regression/SSF/iSSF using glmnet
#'
#' @param f formula
#' @param data data
#'
#' @export
net_clogit <- function(f, data,
                       alpha = 1,
                       type.measure = "deviance",
                       standardize = TRUE,
                       na.action = "na.pass", ...) {

  # NA option
  options(na.action = na.action)

  # extract strata variable and other columns
  wcols <- extract_response_strata(f, other_vars = TRUE)

  # formula with no intercept
  ff <- as.formula(paste0(wcols$case, " ~ -1 + ", wcols$other_vars))
  # explanatory variables
  X <- model.matrix(ff, data)
  # case variable
  p <- data[[wcols$response]]
  # strata
  strata <- data[[wcols$strata]]

  # check data before fitting
  if (anyNA(data)) stop("NA values in data table. Please remove them and rerun.")
  lp <- split(p, strata)
  test <- range(unlist(lapply(lp, mean)))
  if (anyNA(test)) stop("NA values in the response. Please remove them and rerun.")
  if (test[1]==0) stop("There are strata without event. Please remove them and rerun.")
  if (test[1]==1) stop("There are strata without control. Please remove them and rerun.")

  # response variable
  Y <- glmnet::stratifySurv(survival::Surv(rep(1, length(p)), p), strata)

  # fit the model
  fit <- glmnet::glmnet(X, Y, family = "cox",
                        alpha = alpha,
                        type.measure = type.measure,
                        standardize = standardize,
                        ...)

  # return the fitted model
  return(fit)
}

#' Function with similar name
#' @rdname net_clogit
#' @export
net_ssf <- net_clogit

#' Function with similar name
#' @rdname net_clogit
#' @export
net_issf <- net_clogit
