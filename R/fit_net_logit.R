#' Fit logistic regression/RSF with penalized regression in a train-validate-test setup
#'
#' @param ... Options for net_logit and glmnet
#'
#' @rdname fit_net_logit
#' @export
fit_net_logit <- function(f, data,
                          samples, i = 1,
                          metric = c(conditionalBoyce, somersD, AUC)[[1]],
                          out_dir_file = NULL,
                          na.action = "na.pass",
                          ...) {

  # should the sampling be here within the function?

  # get variables
  wcols <- extract_response_strata(f, other_vars = TRUE)

  # filter out NAs and strata with only 1s or 0s
  if(anyNA(data[[wcols$response]])) {
    warning("NAs detected in the response variable. Removing them for model fitting.")
    data <- data[!is.na(data[[wcols$response]]),]
  }

  # separate data for fitting, calibration, and validation
  train_data  <- data[samples$train[[i]], ]
  test_data <- data[samples$test[[i]], ]
  validate_data <- data[samples$validate[[i]], ]

  # perform penalized regression with glmnet
  # use glmnet.cv?
  fit <- net_logit(f, train_data,
                   alpha = 1,
                   type.measure = "deviance",
                   na.action = na.action,
                   ...)

  # get variables
  # f2 <- as.formula(paste0(wcols$response, " ~ -1 + ", wcols$other_vars)) # should we remove the intercept?
  f2 <- as.formula(paste0(wcols$response, " ~ ", wcols$other_vars))

  #----
  # Variable selection step

  # Calibration with conditionalBoyce index (recalled metric here)
  # does not work fro cv.glmnet
  # predict.glmnet(fit, test_data, type = "response")??
  pred_vals <- model.matrix(f2, test_data) %*% coef(fit) # multiple fits?
  # compute conditional Boyce for all lambda parameters
  # here it is just Boyce

  # strat <- rep(1, nrow(test_data))
  # presence <- test_data[[wcols$response]] == 1
  # ecospat::ecospat.boyce(unname(pred_vals[,2])[presence], obs = test_data[[wcols$response]][presence], nclass = 50)
  # conditionalBoyce(data.frame(x = pred_vals[,5], y = test_data[[wcols$response]], strat = rep(1, nrow(test_data))))
  d <- apply(pred_vals, 2, function(x = x, y = y, strat = strat){
    metric(data.frame(x = x, y = y, strat = strat), errors=F)},
    y = test_data[[wcols$response]], strat = rep(1, nrow(test_data)))

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
                                         strat = rep(1, nrow(test_data))))
  results$calibration_score <- max(d)

  #----
  # Validation step
  val_pred_vals <- model.matrix(f2, validate_data) %*% results$coef
  val <- data.frame(x = val_pred_vals,
                    y = validate_data[[wcols$response]],
                    strat = rep(1, nrow(test_data)))
  # val <- split(val, samples$blockH0[sort(match(val$strat, spStrat$id))])
  if(!is.null(samples$blockH0)) {
    val <- split(val, samples$blockH0[match(val$strat, validate_data[[wcols$strata]])])
    results$validation_score <- unlist(lapply(val, metric))
  } else {
    results$validation_score <- metric(val)
  }

  # Validation habitat only
  # no habitat validation score
  results$habitat_validation_score <- NULL
  #plot(results$validation_score, results$habitat_validation_score)

  if (!is.null(out_dir_file)){
    save(results, file = paste0(out_dir_file, "_i", i, ".rda"))
  } else {
    return(results)
  }
}

#' @rdname fit_net_logit
#' @export
fit_net_rsf <- fit_net_logit

#' Fit multiple logistic regression/RSF models with penalized regression in a train-validate-test setup
#'
#' @param ... Options for net_logit and glmnet
#' @param mc.cores Only relevant if `parallel == "mclapply"`. If `parallel == "foreach"`, cores must
#' be assigned before running `fit_multi_net_logit()` using [parallel::makeCluster()] and
#' [doParallel::registerDoParallel()].
#'
#' @rdname fit_multi_net_logit
#' @export
fit_bag_net_logit <- function(f, data,
                                samples,
                                metric = c(conditionalBoyce, somersD, AUC)[[1]],
                                out_dir_file = NULL,
                                na.action = "na.pass",
                                parallel = c(FALSE, "foreach", "mclapply")[1],
                                mc.cores = 2L,
                                ...) {

  # If there is parallel implementation with forach
  if(parallel == "foreach") {
    packs <- c("parallel", "foreach", "doParallel")
    if(!any(packs %in% (base::.packages())))
      warnings(paste0("Parallel fitting of the models using 'foreach' requires the packages ", paste(packs, collapse = ","),
                      " to be loaded and cores to be assigned. Please check it."))
    # check if cores were assigned
    fittedl <- foreach::foreach(i = 1:length(samples$train),
                                .packages = "oneimpact") %dopar% {
                                  fit_net_logit(f = f,
                                                data = data,
                                                samples = samples,
                                                i = i,
                                                metric = metric,
                                                out_dir_file = out_dir_file,
                                                na.action = na.action,
                                                ...)
                                }
  }

  # If there is parallel implementation with forach
  if(parallel == "mclapply") {
    packs <- c("parallel")
    if(!any(packs %in% (base::.packages())))
      warnings(paste0("Parallel fitting of the models using 'mclapply' requires the packages ", paste(packs, collapse = ","),
                      " to be loaded and cores to be assigned. Please check it."))
    # check if cores were assigned
    fitted_list <- parallel::mclapply(1:length(samples$train), mc.cores =  mc.cores, function(i) {
      fit_net_logit(f = f,
                    data = data,
                    samples = samples,
                    i = i,
                    metric = metric,
                    out_dir_file = out_dir_file,
                    na.action = na.action,
                    ...)
    })
  }

  # Common loop if parallel = FALSE
  if(parallel == FALSE) {}
  fitted_list <- list()
  for(i in 1:length(samples$train)) {
    fitted_list[[i]] <- fit_net_logit(f = f,
                                      data = data,
                                      samples = samples,
                                      i = i,
                                      metric = metric,
                                      out_dir_file = out_dir_file,
                                      na.action = na.action,
                                      ...)
  }

  names(fitted_list) <- names(samples$train)
  # define new class?
  fitted_list
}
