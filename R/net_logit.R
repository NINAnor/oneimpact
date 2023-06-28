#' Fits a logistic regression/RSF using glmnet
#'
#' @param f formula
#' @param data data
#' @param alpha Default is L1-regularization (Lasso regression), with `alpha = 1`.
#' L2-regularization (Ridge regression) is done with `alpha = 0`, and elastic-net regression
#' is performed for any `alpha` value between `0` and `1`. For more details, see the
#' [glmnet::glmnet()] documentation.
#' @param na.action Default is `"na.pass"`, i.e. rows with NAs are not automatically
#' removed from the `model.matrix` used for fitting.
#'
#' Check option parallel = TRUE from glmnet.
#'
#' @name net_logit
#' @export
net_logit <- function(f, data,
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
  # response variable
  Y <- data[[wcols$response]]

  # check data before fitting
  if (anyNA(data)) stop("NA values in data table. Please remove them and rerun.")
  if (anyNA(Y)) stop("NA values in the response. Please remove them and rerun.")

  # fit the model
  fit <- glmnet::glmnet(X, Y, family = "binomial",
                        alpha = alpha,
                        type.measure = type.measure,
                        standardize = standardize,
                        ...)

  # return the fitted model
  return(fit)
}

#' Function with similar name
#' @rdname net_logit
#' @export
net_rsf <- net_logit

# Include also net_enm?
