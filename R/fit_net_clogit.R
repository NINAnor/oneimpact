#' Fits a conditional logistic regression/SSF/iSSF using glmnet in a train-validate-test setup
#'
#' By default, [fit_net_clogit()] does not standardize predictor variables. If you want numeric variables
#' to be standardized, you can either use `[bag_git_net_clogit()]` with parameter `standardize = TRUE`
#' or provide an already standardized data set as input.
#'
#' @param f `[formula]` \cr Formula of the model to be fitted, with all possible candidate terms.
#' @param data `[data.frame,tibble]` \cr Complete data set to be analyzed.
#' @param samples `[list]` \cr List of samples with at least three elements: train, test,
#' and validate. Each elements might have several elements, each representing
#' the lines of `data` to be sampled for each resample. Typically, this is computed by
#' the function [oneimpact::create_resamples()].
#' @param method `[character="Lasso"]` \cr The penalized regression method used for fitting
#' each model. Default is `method = "Lasso"`, but it could be `method = "Ridge"` or different
#' flavors of `"AdaptiveLasso"` (see details below).
#' @param metric `[function]{AUC, conditionalBoyce, conditionalSomersD, conditionalAUC}` \cr Function
#' representing the metric to evaluate goodness-of-fit. One of conditionalBoyce (Default),
#' conditionalSomersD, and conditionalAUC. A user-defined function might be provided, with a condition that
#' it must be maximized to find the best fit model.
#' @param kernel_vars `[vector,character=c("step_length", "ta")]` \cr Vector of strings with the names of the variables related
#' to the movement kernel, included in the model (for instance, `"step_length"` and `"turning_angle"`)
#' @param standardize `[logical(1)=TRUE]` \cr Logical flag for predictor variable standardization,
#' prior to fitting the model sequence. The coefficients are always returned on the original scale.
#' Default is standardize=TRUE. If variables are in the same units already, you might not wish to
#' standardize them.
#' @param out_dir_file `[character(1)=NULL]` \cr String with the prefix of the file name (and
#' the folder) where the result of each model will be saved. E.g. if `out_dir_file = "output/test_"`,
#' the models will be saved as RDS files names "test_i1.rds", "test_i2.rds", etc, within the
#' folder "output".
#' @param gamma `[numeric(1)=1]{(0.5, 1, 2)}` \cr Gamma is the exponent for defining the vector of
#' penalty weights when `method = "AdaptiveLasso`. This means that the penalties are defined as
#' `penalty.factor = 1/(coef_ridge^gamma)`, where `coef_ridge` are the coefficients of a Ridge regression.
#' Default is `gamma = 1`, but values of 0.5 or 2 could also be tried, as suggested by the
#' authors (Zou et al 2006).
#' @param ... Options for [oneimpact::net_logit()] and [glmnet::glmnet()].
#'
#' @references Zou, H., 2006. The Adaptive Lasso and Its Oracle Properties. Journal of the American Statistical Association 101, 1418–1429. https://doi.org/10.1198/016214506000000735
#'
#' @name fit_net_clogit
#' @export
fit_net_clogit <- function(f, data,
                           samples, i = 1,
                           kernel_vars = c("step_length", "ta"),
                           metric = c(AUC, conditionalBoyce, conditionalSomersD)[[1]],
                           method = c("Lasso", "Ridge", "AdaptiveLasso",
                                      "DistanceDecay-AdaptiveLasso", "DD-AdaptiveLasso",
                                      "OneZOI-AdaptiveLasso", "OZ-AdaptiveLasso",
                                      "Grouped-AdaptiveLasso", "G-AdaptiveLasso",
                                      "HypothesisDriven-AdaptiveLasso", "HD-AdaptiveLasso",
                                      "ElasticNet")[1],
                           alpha = NULL,
                           penalty.factor = NULL,
                           gamma = 1,
                           standardize = c("internal", FALSE)[1],
                           predictor_table = NULL,
                           function_lasso_decay = c(log, function(x) x/1000)[[1]],
                           factor_hypothesis = 1,
                           factor_grouped_lasso = 1,
                           na.action = "na.pass",
                           out_dir_file = NULL,
                           ...) {

  # parameter checks
  sd_options <- c("internal", "external", FALSE)
  if(!(standardize %in% sd_options))
    stop(paste0("Invalid parameter 'standardize'. It should be one of ", paste(sd_options, collapse = ","), "."))
  method_options <- c("Lasso", "Ridge", "AdaptiveLasso",
                      "DistanceDecay-AdaptiveLasso", "DD-AdaptiveLasso",
                      "OneZOI-AdaptiveLasso", "OZ-AdaptiveLasso",
                      "Grouped-AdaptiveLasso", "G-AdaptiveLasso",
                      "HypothesisDriven-AdaptiveLasso", "HD-AdaptiveLasso",
                      "ElasticNet")
  if(!(grepl(paste(method_options, collapse = "|"), method[1], ignore.case = TRUE)))
    stop(paste0("Invalid parameter 'method'. It should be one of ", paste(method_options, collapse = ","), "."))

  # filter out NAs and strata with only 1s or 0s
  data <- filter_na_strata(f, data)

  # get variables
  wcols <- extract_response_strata(f, covars = TRUE)

  # case
  case <- wcols$response
  # get variable referring to strata
  strat <- wcols$strata

  # relevant columns
  all_vars <- all.vars(f)

  # check columns in data
  if(!all(all_vars %in% names(data)))
    stop(paste0("Not all variables in the formula are present in the data. Please check."))

  # separate data for fitting, calibration, and validation
  if(is.null(samples$sp_strat_id)) {
    # filter by row number
    train_data  <- data[data[[strat]] %in% data[data[[case]] == 1,][[strat]][samples$train[[i]]], all_vars]
    test_data <- data[data[[strat]] %in% data[data[[case]] == 1,][[strat]][samples$test[[i]]], all_vars]
    validate_data <- data[data[[strat]] %in% data[data[[case]] == 1,][[strat]][samples$validate[[i]]], all_vars]
  } else {
    train_data  <- data[data[[strat]] %in% samples$train[[i]], all_vars]
    test_data <- data[data[[strat]] %in% samples$test[[i]], all_vars]
    validate_data <- data[data[[strat]] %in% samples$validate[[i]], all_vars]
  }

  # check NAs
  if(anyNA(train_data)) {
    n_bef <- nrow(train_data)
    # train_data <- filter_na_strata(f, train_data)
    train_data <- filter_na_strata(f, na.omit(train_data))
    nNA <- n_bef - nrow(train_data)
    warning(paste0(nNA, " missing observations were removed from the train set. ", nrow(train_data), " observations were kept."))
  }
  if(anyNA(test_data)) {
    n_bef <- nrow(test_data)
    test_data <- filter_na_strata(f, na.omit(test_data))
    nNA <- n_bef - nrow(test_data)
    warning(paste0(nNA, " missing observations were removed from the test set. ", nrow(test_data), " observations were kept."))
  }
  if(anyNA(validate_data)) {
    n_bef <- nrow(validate_data)
    validate_data <- filter_na_strata(f, na.omit(validate_data))
    nNA <- n_bef - nrow(validate_data)
    warning(paste0(nNA, " missing observations were removed from the validate set. ", nrow(validate_data), " observations were kept."))
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
    if(is.null(predictor_table)) {
      methods_pred_table <- c("Decay-AdaptiveLasso", "DD-AdaptiveLasso",
                              "OneZOI-AdaptiveLasso", "OZ-AdaptiveLasso",
                              "Grouped-AdaptiveLasso", "G-AdaptiveLasso",
                              "HypothesisDriven-AdaptiveLasso", "HD-AdaptiveLasso")
      if(grepl(paste0(methods_pred_table, collapse = "|"), method[1], ignore.case = TRUE)) {
        stop("If 'method' is 'DistanceDecay-AdaptiveLasso', 'OneZOI-AdaptiveLasso', 'Grouped-AdaptiveLasso' or 'HypothesisDriven-AdaptiveLasso', the parameter 'predictor_table' must be provided.")
      }
    }

    # if Decay
    if(grepl("Decay|DD", method[1], ignore.case = TRUE)) {

      # formula
      ff <- as.formula(paste0("~ -1 +", wcols$covars))
      covars <- all.vars(ff)
      # model matrix with data
      M <- stats::model.matrix(ff, data)

      # variables and terms
      terms_order <- attributes(M)$assign
      terms_order <- terms_order[terms_order > 0]
      # vars_formula <- rep(covars, times = unname(table(terms_order)))
      # ZOI and nonZOI variables in the model matrix
      mm_is_zoi <- rep(predictor_table$is_zoi, times = unname(table(terms_order)))
      mm_zoi_radius <- rep(predictor_table$zoi_radius, times = unname(table(terms_order)))
      # cbind(colnames(M), vars_formula, vars_is_zoi, mm_zoi_radius)

      # set penalty factor
      penalty.factor <- ifelse(mm_is_zoi, function_lasso_decay(mm_zoi_radius), 1)
      names(penalty.factor) <- colnames(M)

    } else {

      print("Fitting Ridge...")

      # fit
      ridge_fit <- net_clogit(f, train_data,
                             alpha = 0,
                             type.measure = "deviance",
                             standardize = std,
                             na.action = na.action,
                             ...)
      # get variables
      f2 <- as.formula(paste0(wcols$response, " ~ -1 + ", wcols$covars))
      # calibration
      pred_vals <- model.matrix(f2, test_data) %*% coef(ridge_fit) # multiple fits?
      d <- apply(pred_vals, 2, function(x = x, y = y, strat = strat){
        metric(data.frame(x = x, y = y, strat = strat), errors=F)},
        y = test_data[[wcols$response]], strat = rep(1, nrow(test_data)))
      # coefficients
      coef_ridge <- matrix(coef(ridge_fit)[,which.max(d)]) # coefficients

      #---- prepation to standardize coefs
      if(standardize == "internal") {

        print("Standardizing coefs...")

        # covariates summary
        all_vars <- all.vars(f2)[-1]
        classes <- sapply(data[,all_vars], class)
        # numeric variables
        data_summary_num <- as.data.frame(apply(na.omit(as.matrix(data[,all_vars[classes == "numeric"]])), 2, data_summary))
        # character variables - use mode
        data_summary_ch <- as.data.frame(apply(na.omit(as.matrix(data[,all_vars[classes != "numeric"]])), 2, data_summary_char))
        names(data_summary_ch) <- all_vars[classes != "numeric"]
        dat_summ <- cbind(data_summary_num, data_summary_ch)[order(c(which(classes == "numeric"), which(classes != "numeric")))]

        # info from formula
        ff <- as.formula(paste0("~ -1 +", wcols$covars))
        # covariates
        m_covars <- all.vars(ff, unique = F)
        # are they numeric?
        numeric_covs <- (classes == "numeric")
        repeated <- m_covars[which(duplicated(m_covars))]
        rep_times <- ifelse(names(numeric_covs) %in% repeated, 2, 1) ## CORRECT IF THERE ARE MORE THAN TWO TERMS WITH THE SAME VARIABLE
        numeric_covs <- rep(numeric_covs, times = rep_times)
        # model matrix with data
        M <- stats::model.matrix(ff, data)
        # variables and terms
        terms_order <- attributes(M)$assign
        terms_order <- terms_order[terms_order > 0]
        vars_formula <- rep(m_covars, times = unname(table(terms_order)))
        numeric_vars_order <- rep(numeric_covs, times = unname(table(terms_order)))

        # SDs
        sds <- dat_summ
        sds <- sds[rownames(sds) == "sd", colnames(sds) %in% m_covars]
        sds_all <- sds[match(vars_formula, colnames(sds))]
        # sds_all <- unlist(rep(sds, terms_order)) |>
        #   as.numeric()
        sds_all[numeric_vars_order == FALSE] <- 1
        names(sds_all) <- vars_formula
        sds_all <- unlist(sds_all)

        coef_ridge <- to_std(coef_ridge, sds_all)
      }

      print(method[1])

      # if Adaptive Lasso
      if(tolower(method[1]) == "adaptivelasso") {

        print("Fitting AdaptiveLasso...")

        penalty.factor <- 1/(abs(coef_ridge)**gamma)
        penalty.factor[penalty.factor == Inf] <- 999999999 # If there is any infinite coefficient

      } else {

        # if OneZOI
        if(tolower(method[1]) == "onezoi-adaptivelasso" | tolower(method[1]) == "oz-adaptivelasso") {

          print("Fitting One-ZOI AdaptiveLasso...")

          # prepare from predictor table

          # variables and terms
          terms_order <- attributes(M)$assign
          terms_order <- terms_order[terms_order > 0]
          # vars_formula <- rep(covars, times = unname(table(terms_order)))
          # ZOI and nonZOI variables in the model matrix
          mm_is_zoi <- rep(predictor_table$is_zoi, times = unname(table(terms_order)))
          mm_zoi_radius <- rep(predictor_table$zoi_radius, times = unname(table(terms_order)))
          # cbind(rep(predictor_table$variable, times = unname(table(terms_order))), mm_zoi_radius)
          mm_predictor_vars <- rep(predictor_table$variable, times = unname(table(terms_order)))

          # set penalties
          penalty.factor <- 1/(abs(coef_ridge)**gamma)
          # select only the best
          zoi_terms <- unique(mm_predictor_vars[mm_is_zoi == 1])
          for(i in zoi_terms) {
            vals <- penalty.factor[mm_is_zoi == 1 & mm_predictor_vars == i]
            # keep only the minimum
            vals[vals > min(vals, na.rm = TRUE)] <- Inf
            penalty.factor[mm_is_zoi == 1 & mm_predictor_vars == i] <- vals
          }

          penalty.factor[penalty.factor == Inf] <- 999999999 # If there is any infinite coefficient
          # rownames(penalty.factor) <- mm_predictor_vars

        } else {

          # if grouped
          if(tolower(method[1]) == "grouped-adaptivelasso" | tolower(method[1]) == "g-adaptivelasso") {

            print("Fitting Grouped AdaptiveLasso...")

            # prepare from predictor table

            # variables and terms
            terms_order <- attributes(M)$assign
            terms_order <- terms_order[terms_order > 0]
            # vars_formula <- rep(covars, times = unname(table(terms_order)))
            # ZOI and nonZOI variables in the model matrix
            mm_is_zoi <- rep(predictor_table$is_zoi, times = unname(table(terms_order)))
            mm_zoi_radius <- rep(predictor_table$zoi_radius, times = unname(table(terms_order)))
            mm_predictor_vars <- rep(paste0(predictor_table$variable, predictor_table$cumulative), times = unname(table(terms_order)))

            # set penalties
            penalty.factor <- 1/(abs(coef_ridge)**gamma)

            # select only the best
            zoi_terms <- unique(mm_predictor_vars[mm_is_zoi == 1])
            for(i in zoi_terms) {
              vals <- coef_ridge[mm_is_zoi == 1 & mm_predictor_vars == i]
              # sum relative to the vatiation in the group
              vals <- grouped_func(vals, phi_group = factor_grouped_lasso)**gamma
              penalty.factor[mm_is_zoi == 1 & mm_predictor_vars == i] <- vals
            }

            penalty.factor[penalty.factor == Inf] <- 999999999 # If there is any infinite coefficient

          } else {

            # if grouped
            if(tolower(method[1]) == "hypothesisdriven-adaptivelasso" | tolower(method[1]) == "hd-adaptivelasso") {

              print("Fitting HypothesisDriven AdaptiveLasso...")

              # prepare from predictor table

              # variables and terms
              terms_order <- attributes(M)$assign
              terms_order <- terms_order[terms_order > 0]
              # vars_formula <- rep(covars, times = unname(table(terms_order)))
              # ZOI and nonZOI variables in the model matrix
              mm_is_zoi <- rep(predictor_table$is_zoi, times = unname(table(terms_order)))
              mm_zoi_radius <- rep(predictor_table$zoi_radius, times = unname(table(terms_order)))
              mm_predictor_vars <- rep(predictor_table$variable, times = unname(table(terms_order)))

              # set penalties
              expected_negative <- ifelse(mm_is_zoi == 1, -1, 0)
              penalty.factor <- hypothesis_func(coef_ridge, expected_negative, phi_hyp = factor_hypothesis)**gamma
              # cbind(coef_ridge, expected_negative, penalty.factor)

              # zoi_terms <- unique(mm_predictor_vars[mm_is_zoi == 1])
              # for(i in zoi_terms) {
              #   vals <- penalty.factor[mm_is_zoi == 1 & mm_predictor_vars == i]
              #   # keep only the minimum
              #   vals[vals > min(vals, na.rm = TRUE)] <- Inf
              #   penalty.factor[mm_is_zoi == 1 & mm_predictor_vars == i] <- vals
              # }

              penalty.factor[penalty.factor == Inf] <- 999999999

            } else {

              stop("Please choose a valid method for the function.")
            }
          }
        }
      }
    }
  }

  # perform penalized regression with glmnet
  # use glmnet.cv?
  fit <- net_clogit(f, train_data,
                    alpha = alpha,
                    penalty.factor = penalty.factor,
                    type.measure = "deviance",
                    standardize = std,
                    na.action = na.action,
                    ...)

  # get variables
  f2 <- as.formula(paste0(wcols$response, " ~ -1 + ", wcols$covars))

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
  results$train_score <- metric(data.frame(x = train_pred_vals,
                                           y = train_data[[wcols$response]],
                                           strat = train_data[[wcols$strata]]))
  results$test_score <- max(d)

  #----
  # Validation step
  val_pred_vals <- model.matrix(f2, validate_data) %*% results$coef
  val <- data.frame(x = val_pred_vals,
                    y = validate_data[[wcols$response]],
                    strat = validate_data[[wcols$strata]])
  val<-merge(val, data.frame(strat=samples$sp_strat_id, blockH0=samples$blockH0), by="strat", all.x=T, all.y=F)

  if(!is.null(samples$blockH0)) {

    # data[data$strat %in% validate_data[[wcols$strata]],]$herd |> table()
    val2 <- split(val, val$blockH0)
    # val2 <- split(val, samples$blockH0[match(val$strat, validate_data[[wcols$strata]])])
    # val2 <- split(val, samples$blockH0[val$strat])
    if(length(val2) == 0) {
      if(is.null(samples$sp_strat_id)) {
        val2 <- split(val, samples$blockH0[samples$validate[[i]]])
      } else {
        val2 <- split(val, samples$validate[[i]])
      }
    }
    results$validation_score <- unlist(lapply(val2, metric))

  } else {
    results$validation_score <- metric(val)
  }

  # Validation habitat only
  if(kernel_vars[1] != "") {

    pred_vals_kernel <- kernel_prediction(f, validate_data,
                                          kernel_vars = kernel_vars,
                                          coefs = results$coef[,1])

    pred_vals_habitat <- val_pred_vals - pred_vals_kernel # does it make sense??
    hab <- data.frame(x = pred_vals_habitat,
                      y = validate_data[[wcols$response]],
                      strat = validate_data[[wcols$strata]])
    hab<-merge(hab, data.frame(strat=samples$sp_strat_id, blockH0=samples$blockH0), by="strat", all.x=T, all.y=F)

    if(!is.null(samples$blockH0)) {
      hab2 <- split(hab, hab$blockH0)
      # hab2 <- split(hab, samples$blockH0[match(val$strat, validate_data[[wcols$strata]])])
      # hab2 <- split(hab, samples$blockH0[val$strat])
      if(length(val2) == 0) {
        if(is.null(samples$sp_strat_id)) {
          hab2 <- split(hab, samples$blockH0[samples$validate[[i]]])
        } else {
          hab2 <- split(hab, samples$validate[[i]])
        }
      }
      results$habitat_validation_score <- unlist(lapply(hab2, metric))

    } else {
      results$habitat_validation_score <- metric(hab)
    }
    #plot(results$validation_score, results$habitat_validation_score)


  } else {

    # if there is not movement kernel terms, NULL
    results$habitat_validation_score <- NULL

  }

  # whether to save the results externally
  if (!is.null(out_dir_file)){
    # change if there are more than 999 samples
    names_out <- oneimpact:::pretty_seq(1:999)[i]
    saveRDS(results, file = paste0(out_dir_file, "_", names_out, ".rds"))
  }

  return(results)

}

#' @rdname fit_net_clogit
#' @export
fit_net_ssf <- fit_net_clogit

#' @rdname fit_net_clogit
#' @export
fit_net_issf <- fit_net_clogit

#' Fit a bag of conditional logistic regression/SSF/iSSF models with penalized regression in a train-validate-test setup
#'
#' @param ... Options for net_clogit and glmnet
#' @param mc.cores Only relevant if `parallel == "mclapply"`. If `parallel == "foreach"`, cores must
#' be assigned before running `fit_multi_net_logit()` using [parallel::makeCluster()] and
#' [doParallel::registerDoParallel()].
#' @param subset `[vector]` \cr Vector of samples to be run (e.g. `c(1,2,4,5)` or `3:10`).
#' By default, all the samples in `samples` are run.
#' @param standardize internal = internal glmnet standaridization, i.e. using glmnet with argument standardize = TRUE.
#' This also standardizes dummy variables, but returns the estimated coefficients back to the original scale.
#' This however can cause baises in the estimates because of the bias-variance tradeoff that L1 and L1 regularization
#' methods try to minimize.
#' See more info in https://stackoverflow.com/questions/17887747/how-does-glmnets-standardize-argument-handle-dummy-variables
#' external = glmnet is called with argument standardize = FALSE, but standization is done by the
#' bag_fit_net_logit function. Return coefs in the original scale?? Implement.
#' If FALSE, no standardization of predictors is done.
#'
#' @name bag_fit_net_clogit
#' @export
bag_fit_net_clogit <- function(f, data,
                               samples,
                               subset_samples = 1:length(samples$train),
                               kernel_vars = c("step_length", "ta"),
                               metric = c(conditionalBoyce, conditionalSomersD, conditionalAUC)[[1]],
                               standardize = c("internal", "external", FALSE)[1],
                               method = c("Lasso", "Ridge", "AdaptiveLasso", "DecayAdaptiveLasso", "ElasticNet")[1],
                               alpha = NULL,
                               penalty.factor = NULL,
                               predictor_table = NULL,
                               na.action = "na.pass",
                               out_dir_file = NULL,
                               parallel = c(FALSE, "foreach", "mclapply")[1],
                               mc.cores = 2L,
                               verbose = FALSE,
                               ...) {

  # get variables
  wcols <- extract_response_strata(f, covars = TRUE)

  # First we standardize covariates
  # relevant columns
  all_vars <- all.vars(f)
  all_covars <- all_vars[-1]

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
    fitted_list <- foreach::foreach(i = subset_samples,
                                    .packages = "oneimpact") %dopar% {
                                      try(fit_net_clogit(f = f,
                                                         data = data,
                                                         samples = samples,
                                                         i = i,
                                                         kernel_vars = kernel_vars,
                                                         metric = metric,
                                                         method = method,
                                                         standardize = standardize,
                                                         alpha = alpha,
                                                         penalty.factor = penalty.factor,
                                                         predictor_table = predictor_table,
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
    fitted_list <- parallel::mclapply(subset_samples, function(i) {
      try(fit_net_clogit(f = f,
                         data = data,
                         samples = samples,
                         i = i,
                         kernel_vars = kernel_vars,
                         metric = metric,
                         method = method,
                         standardize = standardize,
                         alpha = alpha,
                         penalty.factor = penalty.factor,
                         predictor_table = predictor_table,
                         na.action = na.action,
                         out_dir_file = out_dir_file,
                         ...))
    }, mc.cores =  mc.cores)
  }

  # Common loop if parallel = FALSE
  if(parallel == FALSE) {
    fitted_list <- list()
    for(i in subset_samples) {
      if(verbose) print(paste0("Fitting sample ", i, "/", length(samples$train), "..."))
      fitted_list[[i]] <- try(fit_net_clogit(f = f,
                                             data = data,
                                             samples = samples,
                                             i = i,
                                             kernel_vars = kernel_vars,
                                             metric = metric,
                                             method = method,
                                             standardize = standardize,
                                             alpha = alpha,
                                             penalty.factor = penalty.factor,
                                             predictor_table = predictor_table,
                                             na.action = na.action,
                                             out_dir_file = out_dir_file,
                                             ...))
    }
  }

  if(length(fitted_list) == results$n)
    names(fitted_list) <- names(samples$train)
  # define new class?
  results$models <- fitted_list

  ## TO DO
  # unstandardize coeffients if standarize = "external"

  # Add info about the covariates - type
  results$numeric_covs <- numeric_covs

  results
}

