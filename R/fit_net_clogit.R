#' Fits a conditional logistic regression/SSF/iSSF using glmnet in a train-validate-test setup
#'
#' @param ... Options for net_clogit and glmnet
#'
#' @name fit_net_clogit
#' @export
fit_net_clogit <- function(f, data,
                           samples, i = 1,
                           kernel_vars = c("step_length", "ta"),
                           metric = c(conditionalBoyce, somersD, AUC)[[1]],
                           out_dir_file = NULL,
                           na.action = "na.pass",
                           ...) {

  # should the sampling be here within the function?

  # filter out NAs and strata with only 1s or 0s
  data <- filter_na_strata(f, data)
  # case
  vars <- extract_response_strata(f)
  case <- vars$response
  # get variable referring to strata
  strat <- vars$strata

  # separate data for fitting, calibration, and validation
  train_data  <- data[data[[strat]] %in% data[data[[case]] == 1,][[strat]][samples$train[[i]]], ]
  test_data <- data[data[[strat]] %in% data[data[[case]] == 1,][[strat]][samples$test[[i]]], ]
  validate_data <- data[data[[strat]] %in% data[data[[case]] == 1,][[strat]][samples$validate[[i]]], ]

  # perform penalized regression with glmnet
  # use glmnet.cv?
  fit <- net_clogit(f, train_data,
                    alpha = 1,
                    type.measure = "deviance",
                    na.action = na.action,
                    ...)

  # get variables
  wcols <- extract_response_strata(f, other_vars = TRUE)
  f2 <- as.formula(paste0(wcols$response, " ~ -1 + ", wcols$other_vars))

  #----
  # Variable selection step

  # Calibration with conditionalBoyce index (recalled metric here)
  # does not work fro cv.glmnet
  # predict.glmnet(fit, test_data, type = "response")??
  pred_vals <- model.matrix(f2, test_data) %*% coef(fit) # multiple fits?
  # compute conditional Boyce for all lambda parameters
  d <- apply(pred_vals, 2, function(x = x, y = y, strat = strat){
    metric(data.frame(x = x, y = y, strat = strat), errors=F)},
    y = test_data[[wcols$response]], strat = test_data[[wcols$strata]])

  # best lambda
  #plot(fit$lambda, d)
  fit$lambda[which.max(d)]

  # initiate results object
  results <- list()
  results$lambda <- fit$lambda[which.max(d)] # lambda
  results$coef <- matrix(coef(fit)[,which.max(d)]) # coefficients
  rownames(results$coef) <- names(coef(fit)[,which.max(d)])
  results$var_names <- names(coef(fit)[,which.max(d)]) # variable names

  # get predicted values based on the training and testing data
  train_pred_vals <- model.matrix(f2, train_data) %*% results$coef
  test_pred_vals <- model.matrix(f2, test_data) %*% results$coef

  # save results
  results$fit_score <- metric(data.frame(x = train_pred_vals,
                                                   y = train_data[[wcols$response]],
                                                   strat = train_data[[wcols$strata]]))
  results$calibration_score <- max(d)

  #----
  # Validation step
  val_pred_vals <- model.matrix(f2, validate_data) %*% results$coef
  val <- data.frame(x = val_pred_vals,
                    y = validate_data[[wcols$response]],
                    strat = validate_data[[wcols$strata]])
  # val <- split(val, samples$blockH0[sort(match(val$strat, spStrat$id))])
  val <- split(val, samples$blockH0[match(val$strat, validate_data[[wcols$strata]])])
  results$validation_score <- unlist(lapply(val, metric))

  # Validation habitat only
  pred_vals_kernel <- kernel_prediction(f, validate_data,
                                        kernel_vars = kernel_vars,
                                        coefs = results$coef[,1])

  pred_vals_habitat <- val_pred_vals - pred_vals_kernel # does it make sense??
  hab <- data.frame(x = pred_vals_habitat,
                    y = validate_data[[wcols$response]],
                    strat = validate_data[[wcols$strata]])
  hab <- split(hab, samples$blockH0[match(hab$strat, validate_data[[wcols$strata]])])
  results$habitat_validation_score <- unlist(lapply(hab, metric))
  #plot(results$validation_score, results$habitat_validation_score)

  if (!is.null(out_dir_file)){
    save(results, file = paste0(out_dir_file, "_i", i, ".rda"))
  } else {
    return(results)
  }
}

#' @rdname fit_net_clogit
#' @export
fit_net_ssf <- fit_net_clogit

#' @rdname fit_net_clogit
#' @export
fit_net_issf <- fit_net_clogit
