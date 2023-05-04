#' Fit issf with spatial stratification and penalized regression
#'
#' @export
spat_strat_fit_issf <- function(f, data, spStrat, sampling_params, i = 0,
           kernel_vars = c("step_length", "ta"), out_dir = NULL) {
    # set points as train, test, and validation
    spStrat$set <- set_validation_sampling(spStrat = spStrat, seed = i,
                                           sampling_params = sampling_params)
    # get only sampled points
    spStrat <- spStrat[!is.na(spStrat$set),]

    # filter out NAs and strata with only 1s or 0s
    data <- filter_na_strata(f, data)
    # get variable referring to strata
    strat <- extract_case_strata(f)$strata
    ###### get number of positions per stratum and filter by that???

    # separate data for fitting, calibration, and validation
    data_fit <- data[data[[strat]] %in% spStrat$id[spStrat$set == "fitting"],]
    data_cal <- data[data[[strat]] %in% spStrat$id[spStrat$set == "calibration"],]
    data_val <- data[data[[strat]] %in% spStrat$id[spStrat$set == "validation"],]

    # perform penalized regression with glmnet
    # use glmnet.cv?
    fit <- net_issf(f, data_fit)

    # get variables
    wcols <- extract_case_strata(f, other_vars = TRUE)
    f2 <- as.formula(paste0(wcols$case, " ~ -1 + ", wcols$other_vars))

    # Calibration with conditionalBoyce index
    # does not work fro cv.glmnet
    predVals <- model.matrix(f2, data_cal) %*% coef(fit) # multiple fits?
    #predict.glmnet(fit, data_cal, type = "response")??
    # compute conditional Boyce for all lambda parameters
    d <- apply(predVals, 2, function(x = x, y = y, strat = strat){
      conditionalBoyce(data.frame(x = x, y = y, strat = strat), errors=F)},
      y = data_cal[[extract_case_strata(f)$case]], strat = data_cal[[extract_case_strata(f)$strata]])

    #plot(fit$lambda, d)
    fit$lambda[which.max(d)]

    # get predicted values based on the training and testing data
    fitVals <- model.matrix(f2, data_fit) %*% matrix(coef(fit)[,which.max(d)])
    predVals <- model.matrix(f2, data_cal) %*% matrix(coef(fit)[,which.max(d)])

    # save results
    results <- list()
    results$coef <- matrix(coef(fit)[,which.max(d)]) # coefficients
    results$var_names <- names(coef(fit)[,which.max(d)]) # variable names
    results$lambda <- fit$lambda[which.max(d)] # lambda
    results$fit_score <- conditionalBoyce(data.frame(x = fitVals,
                                                     y = data_fit[[extract_case_strata(f)$case]],
                                                     strat = data_fit[[extract_case_strata(f)$strata]]))
    results$calibration_score <- max(d)

    # Validation
    predVals <- model.matrix(f2, data_val) %*% results$coef
    test <- data.frame(x = predVals,
                       y = data_val[[extract_case_strata(f)$case]],
                       strat = data_val[[extract_case_strata(f)$strata]])
    test <- split(test, spStrat$blockH0[match(test$strat, spStrat$id)])
    results$validation_score <- unlist(lapply(test, conditionalBoyce))

    # Validation habitat only
    predVals_kernel <- kernel_prediction(f, data_val,
                                         kernel_vars = kernel_vars,
                                         coefs=coef(fit)[,which.max(d)])

    predVals_habitat <- predVals - predVals_kernel # does it make sense??
    test <- data.frame(x = predVals_habitat, y = data_val[[extract_case_strata(f)$case]],
                       strat = data_val[[extract_case_strata(f)$strata]])
    test <- split(test, spStrat$blockH0[match(test$strat, spStrat$id)])
    results$habitat_validation_score <- unlist(lapply(test, conditionalBoyce))
    #plot(results$validation_score, results$habitat_validation_score)

    if (!is.null(out_dir)){
      save(results, file=paste0(out_dir, "spat_strat_issf_i", i, ".rda"))
    }else{
      return(results)
    }
  }
