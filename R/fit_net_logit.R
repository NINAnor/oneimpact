#' Fit logistic regression/RSF with penalized regression in a train-validate-test setup
#'
#' By default, [fit_net_logit()] does not standardize predictor variables. If you want numeric variables
#' to be standardized, you can either use `[bag_git_net_logit()]` with parameter `standardize = TRUE`
#' or provide an already standardized data set as input.
#'
#' @param f `[formula]` \cr Formula of the model to be fitted, with all possible candidate terms.
#' @param data `[data.frame,tibble]` \cr Complete data set to be analyzed.
#' @param samples `[list]` \cr List of samples with at least three elements: train, test,
#' and validate. Each elements might have several elements, each representing
#' the lines of `data` to be sampled for each resample. Typically, this is computed by
#' the function [oneimpact::create_resamples()].
#' @param metric `[function]{conditionalBoyce, SomersD, AUC, proc_AUC}` \cr Function
#' representing the metric to evaluate goodness-of-fit. One of conditionalBoyce (Default),
#' somersD, AUC, and proc_AUC. A user-defined function might be provided, with a condition that
#' it must be maximized to find the best fit model.
#' @param standardize `[logical(1)=TRUE]` \cr Logical flag for predictor variable standardization,
#' prior to fitting the model sequence. The coefficients are always returned on the original scale.
#' Default is standardize=TRUE. If variables are in the same units already, you might not wish to
#' standardize them.
#' @param out_dir_file `[character(1)=NULL]` \cr String with the prefix of the file name (and
#' the folder) where the result of each model will be saved. E.g. if `out_dir_file = "output/test_"`,
#' the models will be saved as RDS files names "test_i1.rds", "test_i2.rds", etc, within the
#' folder "output".
#' @param ... Options for net_logit and glmnet
#'
#' @name fit_net_logit
#' @export
fit_net_logit <- function(f, data,
                          samples, i = 1,
                          metric = c(conditionalBoyce, somersD, AUC, proc_AUC)[[1]],
                          method = c("Lasso", "Rigdge", "AdaptiveLasso", "DecayAdaptiveLasso", "ElasticNet")[1],
                          alpha = NULL,
                          penalty.factor = NULL,
                          standardize = c("internal", FALSE)[1],
                          predictor_grid = NULL,
                          lasso_decay_type = c(log, function(x) x/1000)[[1]],
                          na.action = "na.pass",
                          out_dir_file = NULL,
                          ...) {

  # parameter checks
  sd_options <- c("internal", "external", FALSE)
  if(!(standardize %in% sd_options))
    stop(paste0("Invalid parameter 'standardize'. It should be one of ", paste(sd_options, collapse = ","), "."))
  method_options <- c("Lasso", "Ridge", "AdaptiveLasso", "DecayAdaptiveLasso", "ElasticNet")
  if(!(grepl(paste(method_options, collapse = "|"), method[1], ignore.case = TRUE)))
    stop(paste0("Invalid parameter 'method'. It should be one of ", paste(method_options, collapse = ","), "."))

  # get variables
  wcols <- extract_response_strata(f, other_vars = TRUE)

  # filter out NAs
  if(anyNA(data[[wcols$response]])) {
    warning("NAs detected in the response variable. Removing them for model fitting.")
    data <- data[!is.na(data[[wcols$response]]),]
  }

  # separate data for fitting, calibration, and validation
  train_data  <- data[samples$train[[i]], ]
  test_data <- data[samples$test[[i]], ]
  validate_data <- data[samples$validate[[i]], ]

  # check NAs
  if(anyNA(train_data)) {
    train_data <- na.omit(train_data)
    nNA <- length(na.action(train_data))
    warning(paste0(nNA, " missing observations were removed from the test set. ", nrow(train_data), " observations were kept."))
  }
  if(anyNA(test_data)) {
    test_data <- na.omit(test_data)
    nNA <- length(na.action(test_data))
    warning(paste0(nNA, " missing observations were removed from the test set. ", nrow(test_data), " observations were kept."))
  }
  if(anyNA(validate_data)) {
    validate_data <- na.omit(validate_data)
    nNA <- length(na.action(validate_data))
    warning(paste0(nNA, " missing observations were removed from the test set. ", nrow(validate_data), " observations were kept."))
  }

  # set standardize parameter to be used in glmnet call
  if(standardize != "internal") {
    std <- FALSE
  } else {
    std <- TRUE
  }

  # set method - alpha
  if(is.null(alpha)) {
    # If Lasso
    if(grepl("Lasso", method[1], ignore.case = TRUE)) {
      alpha <- 1
    } else {
      # If Ridge
      if(grepl("Ridge", method[1], ignore.case = TRUE)) {
        alpha <- 0
      }
    }
  }

  # set method - penalties
  if(is.null(penalty.factor)) {

    # check
    # variable grid to define penalties
    if(is.null(predictor_grid)) {
      if(grepl("DecayAdaptiveLasso", method[1], ignore.case = TRUE)) {
        stop("If 'method' is 'DecayAdaptiveLasso', the parameter 'predictor_grid' must be provided.")
      }
    }

    # if Decay
    if(grepl("Decay", method[1], ignore.case = TRUE)) {

      # formula
      ff <- as.formula(paste0("~ -1 +", wcols$other_vars))
      covars <- all.vars(ff)
      # model matrix with data
      M <- stats::model.matrix(ff, data)

      # variables and terms
      terms_order <- attributes(M)$assign
      terms_order <- terms_order[terms_order > 0]
      # vars_formula <- rep(covars, times = unname(table(terms_order)))
      # ZOI and nonZOI variables in the model matrix
      mm_is_zoi <- rep(predictor_grid$is_zoi, times = unname(table(terms_order)))
      mm_zoi_radius <- rep(predictor_grid$zoi_radius, times = unname(table(terms_order)))
      # cbind(colnames(M), vars_formula, vars_is_zoi, mm_zoi_radius)

      # set penalty factor
      penalty.factor <- ifelse(mm_is_zoi, lasso_decay_type(mm_zoi_radius), 1)
      names(penalty.factor) <- colnames(M)

    } else {
      if(tolower(method[1]) == "adaptivelasso") {

        # fit
        ridge_fit <- net_logit(f, train_data,
                               alpha = 0,
                               type.measure = "deviance",
                               standardize = std,
                               na.action = na.action,
                               ...)
        # get variables
        f2 <- as.formula(paste0(wcols$response, " ~ -1 + ", wcols$other_vars))
        # calibration
        pred_vals <- model.matrix(f2, test_data) %*% coef(ridge_fit)[-1,] # multiple fits?
        d <- apply(pred_vals, 2, function(x = x, y = y, strat = strat){
          metric(data.frame(x = x, y = y, strat = strat), errors=F)},
          y = test_data[[wcols$response]], strat = rep(1, nrow(test_data)))
        coef_weights <- matrix(coef(ridge_fit)[-1,which.max(d)]) # coefficients

        penalty.factor <- 1/coef_weights
      }
    }
  }

  # perform penalized regression with glmnet
  # use glmnet.cv?
  fit <- net_logit(f, train_data,
                   alpha = alpha,
                   penalty.factor = penalty.factor,
                   type.measure = "deviance",
                   standardize = std,
                   na.action = na.action,
                   ...)

  # get variables
  f2 <- as.formula(paste0(wcols$response, " ~ -1 + ", wcols$other_vars)) # should we remove the intercept?
  #f2 <- as.formula(paste0(wcols$response, " ~ ", wcols$other_vars))

  #----
  # Variable selection step

  # Calibration with conditionalBoyce index (recalled metric here)
  # does not work fro cv.glmnet
  # predict.glmnet(fit, test_data, type = "response")??
  pred_vals <- model.matrix(f2, test_data) %*% coef(fit)[-1,] # multiple fits?
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

  # initialize results
  results <- list()
  results$lambda <- fit$lambda[which.max(d)] # lambda
  results$coef <- matrix(coef(fit)[-1,which.max(d)]) # coefficients
  rownames(results$coef) <- names(coef(fit)[-1,which.max(d)])
  results$var_names <- names(coef(fit)[-1,which.max(d)]) # variable names

  # get predicted values based on the training and testing data
  train_pred_vals <- model.matrix(f2, train_data) %*% results$coef
  test_pred_vals <- model.matrix(f2, test_data) %*% results$coef

  # save results
  results$train_score <- metric(data.frame(x = train_pred_vals,
                                           y = train_data[[wcols$response]],
                                           strat = rep(1, nrow(train_data))))
  results$test_score <- max(d)

  #----
  # Validation step
  val_pred_vals <- model.matrix(f2, validate_data) %*% results$coef
  val <- data.frame(x = val_pred_vals,
                    y = validate_data[[wcols$response]],
                    strat = rep(1, nrow(validate_data)))
  # val <- split(val, samples$blockH0[sort(match(val$strat, spStrat$id))])
  if(!is.null(samples$blockH0)) {

    val2 <- split(val, samples$blockH0[match(val$strat, validate_data[[wcols$strata]])])
    if(length(val2) == 0) {
      val2 <- split(val, samples$blockH0[samples$validate[[i]]])
    }
    results$validation_score <- unlist(lapply(val2, metric))

  } else {
    results$validation_score <- metric(val)
  }

  # Validation habitat only
  # no habitat validation score
  results$habitat_validation_score <- NULL
  #plot(results$validation_score, results$habitat_validation_score)

  # whether to save the results externally
  if (!is.null(out_dir_file)){
    saveRDS(results, file = paste0(out_dir_file, "_i", i, ".rds"))
  } else {
    return(results)
  }
}

#' @rdname fit_net_logit
#' @export
fit_net_rsf <- fit_net_logit

#' Fit a bag of logistic regression/RSF models with penalized regression in a train-validate-test setup
#'
#' @param ... Options for net_logit and glmnet
#' @param mc.cores Only relevant if `parallel == "mclapply"`. If `parallel == "foreach"`, cores must
#' be assigned before running `fit_multi_net_logit()` using [parallel::makeCluster()] and
#' [doParallel::registerDoParallel()].
#' @param standardize internal = internal glmnet standaridization, i.e. using glmnet with argument standardize = TRUE.
#' This also standardizes dummy variables, but returns the estimated coefficients back to the original scale.
#' This however can cause baises in the estimates because of the bias-variance tradeoff that L1 and L1 regularization
#' methods try to minimize.
#' See more info in https://stackoverflow.com/questions/17887747/how-does-glmnets-standardize-argument-handle-dummy-variables
#' external = glmnet is called with argument standardize = FALSE, but standization is done by the
#' bag_fit_net_logit function. Return coefs in the original scale?? Implement.
#' If FALSE, no standardization of predictors is done.
#'
#' @name bag_fit_net_logit
#' @export
bag_fit_net_logit <- function(f, data,
                              samples,
                              metric = c(conditionalBoyce, somersD, AUC, proc_AUC)[[1]],
                              standardize = c("internal", "external", FALSE)[1],
                              method = c("Lasso", "Rigdge", "AdaptiveLasso", "DecayAdaptiveLasso", "ElasticNet")[1],
                              alpha = NULL,
                              penalty.factor = NULL,
                              predictor_grid = NULL,
                              na.action = "na.pass",
                              out_dir_file = NULL,
                              parallel = c(FALSE, "foreach", "mclapply")[1],
                              mc.cores = 2L,
                              verbose = FALSE,
                              ...) {

  # get variables
  wcols <- extract_response_strata(f, other_vars = TRUE)

  # First we standardize covariates
  # relevant columns
  all_vars <- all.vars(f)
  all_covars <- all_vars[grep(wcols$response, all_vars, invert = TRUE)]

  # get predictors
  data_covs <- data[, all_covars]
  # select numeric predictors to be standardized
  numeric_covs <- (sapply(data_covs, class) == "numeric")
  # standardize
  if(standardize == "external") {
    data_covs_num <- data_covs[, numeric_covs]
    # standardize
    data_covs_num_std <- lapply(1:ncol(data_covs_num), function(i) scale(data_covs_num[,i]))
    # register mean and sd
    covs_mean_sd <- data.frame(do.call("rbind",lapply(1:length(data_covs_num_std), function(i)
      sapply(c("scaled:center", "scaled:scale"), function(m) attr(data_covs_num_std[[i]], m)))))
    rownames(covs_mean_sd) <- colnames(data_covs_num)
    colnames(covs_mean_sd) <- c("mean", "sd")
    # merge standardized predictors with non numeric predictors
    data_covs_std <- cbind(data_covs[, !numeric_covs], data.frame(do.call("cbind", data_covs_num_std)))
    data_covs_std <- data_covs_std[,order(c(which(!numeric_covs), which(numeric_covs)))]
    colnames(data_covs_std) <- colnames(data_covs)
    data <- cbind(data[wcols$response], data_covs_std)
  } else {
    data <- data[, all_vars]
  }

  # initiate results object
  results <- list()
  results$n <- length(samples$train)
  results$formula <- f
  results$method <- method
  results$metric <- metric
  # standarized means and sd
  if(standardize == "external") {
    results$covariate_mean_sd <- covs_mean_sd
  } else {
    results$covariate_mean_sd <- NULL
  }

  # If there is parallel implementation with forach
  if(parallel == "foreach") {
    packs <- c("parallel", "foreach", "doParallel")
    if(!any(packs %in% (base::.packages())))
      warnings(paste0("Parallel fitting of the models using 'foreach' requires the packages ", paste(packs, collapse = ","),
                      " to be loaded and cores to be assigned. Please check it."))
    # check if cores were assigned
    fittedl <- foreach::foreach(i = 1:length(samples$train),
                                .packages = "oneimpact") %dopar% {
                                  try(fit_net_logit(f = f,
                                                    data = data,
                                                    samples = samples,
                                                    i = i,
                                                    metric = metric,
                                                    method = method,
                                                    standardize = standardize,
                                                    alpha = alpha,
                                                    penalty.factor = penalty.factor,
                                                    predictor_grid = predictor_grid,
                                                    na.action = na.action,
                                                    out_dir_file = out_dir_file,
                                                    ...))
                                }
  }

  # If there is parallel implementation with forach
  if(parallel == "mclapply") {
    packs <- c("parallel")
    if(!any(packs %in% (base::.packages())))
      warnings(paste0("Parallel fitting of the models using 'mclapply' requires the packages ", paste(packs, collapse = ","),
                      " to be loaded and cores to be assigned. Please check it."))
    # check if cores were assigned
    fitted_list <- parallel::mclapply(1:length(samples$train), function(i) {
      try(fit_net_logit(f = f,
                        data = data,
                        samples = samples,
                        i = i,
                        metric = metric,
                        method = method,
                        standardize = standardize,
                        alpha = alpha,
                        penalty.factor = penalty.factor,
                        predictor_grid = predictor_grid,
                        na.action = na.action,
                        out_dir_file = out_dir_file,
                        ...))
    }, mc.cores =  mc.cores)
  }

  # Common loop if parallel = FALSE
  if(parallel == FALSE) {
    fitted_list <- list()
    for(i in 1:length(samples$train)) {
      if(verbose) print(paste0("Fitting sample ", i, "/", length(samples$train), "..."))
      fitted_list[[i]] <- try(fit_net_logit(f = f,
                                            data = data,
                                            samples = samples,
                                            i = i,
                                            metric = metric,
                                            method = method,
                                            standardize = standardize,
                                            alpha = alpha,
                                            penalty.factor = penalty.factor,
                                            predictor_grid = predictor_grid,
                                            na.action = na.action,
                                            out_dir_file = out_dir_file,
                                            ...))
    }
  }

  names(fitted_list) <- names(samples$train)
  # define new class?
  results$models <- fitted_list

  # Add info about the covariates - type
  results$numeric_covs <- numeric_covs

  ##############
  # if standardize == "external", we should unstandardize the coefs here!
  # implement that later

  results
}
