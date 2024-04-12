#' Fits a logistic regression/RSF using glmnet
#'
#' @param f `[formula]` \cr Formula of the model to be fitted, with all possible candidate terms.
#' @param data `[data.frame,tibble]` \cr Complete data set to be analyzed.#' @param alpha Default is L1-regularization (Lasso regression), with `alpha = 1`.
#' L2-regularization (Ridge regression) is done with `alpha = 0`, and elastic-net regression
#' is performed for any `alpha` value between `0` and `1`. For more details, see the
#' [glmnet::glmnet()] documentation. For Adaptive and Decay Adaptive Lasso, keep `alpha = 1`.
#' @param penalty.factor `[numeric,vector=NULL]` \cr Vector of penalty factors to be used for Adaptive Lasso
#' fitting. The vector might have the same length as the the number of columns given by the model matrix,
#' `model.matrix(f, data)`. Default is `NULL`, in case the same penalty is applied to all variables.
#' @param type.measure `[character(1)="deviance"]` \cr Type of measure to evaluate the model internally
#' in [glmnet::glmnet()]. For logistic and conditional logistic regression, it is by default `"deviance"`.
#' @param standardize `[logical(1)=TRUE]` \cr Whether the call to `glmnet` should include
#' internal standardization of variables or not. Default is TRUE. It should be set to `FALSE`
#' if the variables are already standardized.
#' @param na.action `[character(1)="na.pass"]` \cr Default is `"na.pass"`, i.e. rows with NAs are not automatically
#' removed from the `model.matrix` used for fitting.
#' @param func `[character(1)="glmnet"]{"glmnet", "cv.glmnet"}` \cr The function to be used for
#' fitting. Default is [glmnet::glmnet()]. The second option is [glmnet::cv.glmnet()] which
#' already performs the cross-validation and might include the variable selection/callibration
#' within.
#'
#' Check option parallel = TRUE from glmnet.
#'
#' @name net_logit
#' @export
net_logit <- function(f, data,
                      alpha = 1,
                      penalty.factor = NULL,
                      type.measure = "deviance",
                      standardize = TRUE,
                      na.action = "na.pass",
                      func = c("glmnet", "cv.glmnet")[1],
                      ...) {

  # NA option
  options(na.action = na.action)

  # extract strata variable and other columns
  wcols <- extract_response_strata(f, covars = TRUE)

  # formula with no intercept
  ff <- as.formula(paste0(wcols$reponse, " ~ -1 + ", wcols$covars))
  # explanatory variables
  X <- model.matrix(ff, data)
  # response variable
  Y <- data[[wcols$response]]

  # check data before fitting
  if (anyNA(data)) stop("NA values in data table. Please remove them and rerun.")
  if (anyNA(Y)) stop("NA values in the response. Please remove them and rerun.")

  # set penalty factors if NULL
  if(is.null(penalty.factor))
    penalty.factor <- rep(1, ncol(X))

  # fit the model
  if(func == "glmnet") {
    fit <- glmnet::glmnet(X, Y, family = "binomial",
                          alpha = alpha,
                          penalty.factor = penalty.factor,
                          type.measure = type.measure,
                          standardize = standardize,
                          ...)
  } else {
    fit <- glmnet::cv.glmnet(X, Y, family = "binomial",
                             alpha = alpha,
                             penalty.factor = penalty.factor,
                             type.measure = type.measure,
                             standardize = standardize,
                             ...)
  }

  # return the fitted model
  return(fit)
}

#' Function with similar name
#' @rdname net_logit
#' @export
net_rsf <- net_logit

# Include also net_enm?
