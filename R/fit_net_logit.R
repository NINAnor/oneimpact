#' Fit logistic regression/RSF with penalized regression using glmnet in a train-validate-test setup
#'
#' By default, [fit_net_logit()] does not standardize predictor variables. If you want numeric variables
#' to be standardized, you can either use `[bag_fit_net_logit()]` with parameter `standardize = TRUE`
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
#' @param metric `[function,character]{AUC, conditionalBoyce, conditionalSomersD, conditionalAUC}` \cr Function
#' representing the metric to evaluate goodness-of-fit. One of AUC (Default), conditionalBoyce,
#' conditionalSomersD, and conditionalAUC. A user-defined function might be provided, with a condition that
#' it must be maximized to find the best fit model. It can also be a character, in case it should be one
#' of the following: `c("AUC", "conditionalAUC", "conditionalBoyce", "conditionalSomersD")`.
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
#' @param replace_missing_NA `[logical(1)=TRUE]` \cr If `TRUE` (default), any variables missing from the data
#' (i.e. with variance zero) are removed from the formula for the model fitting procedure, and a `NA` is set
#' as its coefficient in the output. If `FALSE`, the function raises an error if there are variables with
#' variance zero in the formula.
#' @param ... Options for [oneimpact::net_logit()] and [glmnet::glmnet()].
#'
#' @references Zou, H., 2006. The Adaptive Lasso and Its Oracle Properties. Journal of the American Statistical Association 101, 1418–1429. https://doi.org/10.1198/016214506000000735
#'
#' @name fit_net_functions
#' @export
fit_net_logit <- function(f, data,
                          samples, i = 1,
                          metric = c("AUC")[1],
                          metrics_evaluate = c("AUC"),
                          method = c("Lasso", "Ridge", "AdaptiveLasso",
                                     "DistanceDecayLasso", "DDLasso",
                                     "TruncatedLasso",
                                     "OneZOI-AdaptiveLasso", "OZ-AdaptiveLasso",
                                     "Grouped-AdaptiveLasso", "G-AdaptiveLasso",
                                     "HypothesisDriven-AdaptiveLasso", "HD-AdaptiveLasso",
                                     "ElasticNet")[1],
                          alpha = NULL,
                          penalty.factor = NULL,
                          gamma = 1,
                          standardize = c("internal", "external", FALSE)[1],
                          predictor_table = NULL,
                          function_lasso_decay = c(log, function(x) x/1000)[[1]],
                          value_lasso_decay = 1,
                          function_hypothesis = c(exp)[[1]],
                          expected_sign_hypothesis = -1,
                          factor_grouped_lasso = 1,
                          replace_missing_NA = TRUE,
                          na.action = "na.pass",
                          out_dir_file = NULL,
                          verbose = FALSE,
                          ...) {

  #-----------------------------
  # record initial parameters
  parms <- list(call = "fit_net_logit",
                f = f,
                samples = samples,
                i = i,
                metric = metric,
                method = method,
                alpha = alpha,
                penalty.factor = penalty.factor,
                gamma = gamma,
                standardize = standardize,
                predictor_table = predictor_table,
                function_lasso_decay = function_lasso_decay,
                value_lasso_decay = value_lasso_decay,
                function_hypothesis = function_hypothesis,
                expected_sign_hypothesis = expected_sign_hypothesis,
                factor_grouped_lasso = factor_grouped_lasso,
                replace_missing_NA = replace_missing_NA)

  #-----------------------------
  # parameter checks
  # standardize
  sd_options <- c("internal", "external", FALSE)
  if(!(standardize %in% sd_options))
    stop(paste0("Invalid parameter 'standardize'. It should be one of ", paste(sd_options, collapse = ","), "."))
  # method
  method_options <- c("Lasso", "Ridge", "AdaptiveLasso",
                      "DistanceDecayLasso", "DDLasso",
                      "TruncatedLasso",
                      "OneZOI-AdaptiveLasso", "OZ-AdaptiveLasso",
                      "Grouped-AdaptiveLasso", "G-AdaptiveLasso",
                      "HypothesisDriven-AdaptiveLasso", "HD-AdaptiveLasso",
                      "ElasticNet")
  if(!(grepl(paste(method_options, collapse = "|"), method[1], ignore.case = TRUE)))
    stop(paste0("Invalid parameter 'method'. It should be one of ", paste(method_options, collapse = ","), "."))
  # metric
  metric_options <- c("coxnet.deviance", "Cindex", "AUC", "conditionalAUC", "conditionalBoyce", "conditionalSomersD")
  if(is.character(metric)) {
    if(!(metric %in% metric_options)) {
      stop(paste0("Invalid parameter 'metric'. It should be one of ", paste(metric_options, collapse = ","), " or a function."))
    } else {
      metric_fun <- getFromNamespace(metric, ns = "oneimpact")
    }
  } else {
    metric_fun <- metric
  }

  #-----------------------------
  # prepare data

  # get variables
  wcols <- extract_response_strata(f, covars = TRUE)

  # case
  case <- wcols$response

  # filter out NAs in the response variable
  if(anyNA(data[[wcols$response]])) {
    warning("NAs detected in the response variable. Removing them for model fitting.")
    data <- data[!is.na(data[[case]]),]
  }

  # relevant columns
  all_vars <- all.vars(f)

  # check columns in data
  if(!all(all_vars %in% names(data)))
    stop(paste0("Not all variables in the formula are present in the data. Please check."))

  #---
  # Need to standardize?

  # First we standardize covariates, if needed
  # in any case, we keep their mean and sd which might be useful
  all_covars <- all_vars[-1]
  # get predictors
  data_covs <- data[, all_covars]
  # select numeric predictors to be standardized
  numeric_covs <- (sapply(data_covs, is.numeric) == TRUE)

  # if(standardize != FALSE) {
  # now we compute the means and sds in all cases, to check
  # for missing variables/infrastructures

  # get standardization parms
  data_covs_num <- data_covs[, numeric_covs]
  # standardize
  data_covs_num_std <- lapply(1:ncol(data_covs_num), function(i) scale(data_covs_num[,i]))
  # register mean and sd
  covs_mean_sd <- data.frame(do.call("rbind",lapply(1:length(data_covs_num_std), function(i)
    sapply(c("scaled:center", "scaled:scale"), function(m) attr(data_covs_num_std[[i]], m)))))
  rownames(covs_mean_sd) <- colnames(data_covs_num)
  colnames(covs_mean_sd) <- c("mean", "sd")

  # check if there are NAs (no variation in a variables)
  any_missing_vars <- FALSE
  if(any(ind <- which(covs_mean_sd$sd == 0))) {
    # stop("There are covariates with variance equals zero; please check.")
    warning(paste0("Some of the covariates have variance zero: ",
                   paste0(rownames(covs_mean_sd)[ind], collapse = ",")))

    any_missing_vars <- TRUE

    # check if we remove missing variables with NA
    if(replace_missing_NA) {

      # backup original formula
      f_original <- f
      # remove from formula
      vv <- 1
      for(vv in seq_along(ind)) {
        f <- update(f, paste0("~ . - ", rownames(covs_mean_sd)[ind[vv]]))
      }
      warning(paste0("The formula was updated excluding these variable(s): ",
                     paste0(rownames(covs_mean_sd)[ind], collapse = ","), ". ",
                     "The model will be fitted with this updated formula."))
      # parms$f <- f
      # update covariates
      wcols <- extract_response_strata(f, covars = TRUE)

    } else {
      stop("Computation interrupted. Please check the covariates.")
    }

  }

  # }

  # register/use standardized covariates
  if(standardize == "external") {
    # merge standardized predictors with non numeric predictors
    data_covs_std <- cbind(data_covs[, !numeric_covs], data.frame(do.call("cbind", data_covs_num_std)))
    data_covs_std <- data_covs_std[,order(c(which(!numeric_covs), which(numeric_covs)))]
    colnames(data_covs_std) <- colnames(data_covs)
    data <- cbind(data[wcols$response], data_covs_std)
  } else {
    # or just use the original data
    data <- data[, all_vars]
  }

  # separate data for fitting, calibration, and validation
  train_data  <- data[samples$train[[i]], all_vars]
  test_data <- data[samples$test[[i]], all_vars]
  validate_data <- data[samples$validate[[i]], all_vars]

  # check NAs
  if(anyNA(train_data)) {
    train_data <- na.omit(train_data)
    nNA <- length(na.action(train_data))
    warning(paste0(nNA, " missing observations were removed from the train set. ", nrow(train_data), " observations were kept."))
  }
  if(anyNA(test_data)) {
    test_data <- na.omit(test_data)
    nNA <- length(na.action(test_data))
    warning(paste0(nNA, " missing observations were removed from the test set. ", nrow(test_data), " observations were kept."))
  }
  if(anyNA(validate_data)) {
    validate_data <- na.omit(validate_data)
    nNA <- length(na.action(validate_data))
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
      } else {
        if(grepl("ElasticNet", method[1], ignore.case = TRUE)) {
          alpha <- 0.5
        }
      }
    }
  }

  # set method - penalties
  if(is.null(penalty.factor)) {

    # check
    # variable grid to define penalties
    if(is.null(predictor_table)) {
      methods_pred_table <- c("DistanceDecayLasso", "DDLasso",
                              "OneZOI-AdaptiveLasso", "OZ-AdaptiveLasso",
                              "Grouped-AdaptiveLasso", "G-AdaptiveLasso",
                              "HypothesisDriven-AdaptiveLasso", "HD-AdaptiveLasso")
      if(grepl(paste0(methods_pred_table, collapse = "|"), method[1], ignore.case = TRUE)) {
        stop("If 'method' is 'DistanceDecayLasso', 'OneZOI-AdaptiveLasso', 'Grouped-AdaptiveLasso' or 'HypothesisDriven-AdaptiveLasso', the parameter 'predictor_table' must be provided.")
      }
    }

    # formula
    ff <- as.formula(paste0("~ -1 +", wcols$covars))
    covars <- all.vars(ff)
    # model matrix with data
    M <- stats::model.matrix(ff, data)

    # if Decay - without the Adaptive part, just Lasso
    if(grepl("Decay|DD", method[1], ignore.case = TRUE)) {

      # variables and terms
      terms_order <- attributes(M)$assign
      terms_order <- terms_order[terms_order > 0]
      # vars_formula <- rep(covars, times = unname(table(terms_order)))
      # ZOI and nonZOI variables in the model matrix
      mm_is_zoi <- rep(predictor_table$is_zoi, times = unname(table(terms_order)))
      mm_zoi_radius <- rep(predictor_table$zoi_radius, times = unname(table(terms_order)))
      # cbind(colnames(M), vars_formula, vars_is_zoi, mm_zoi_radius)

      # set penalty factor
      penalty.factor <- ifelse(mm_is_zoi, function_lasso_decay(mm_zoi_radius), value_lasso_decay)
      names(penalty.factor) <- colnames(M)

    } else {

      if(tolower(method[1]) != "lasso" & tolower(method[1]) != "elasticnet") {

        if(verbose) print("Fitting Ridge...")

        # fit
        ridge_fit <- net_logit(f, train_data,
                               alpha = 0,
                               type.measure = "deviance",
                               standardize = std,
                               na.action = na.action,
                               ...)

        if(tolower(method[1]) != "ridge") {

          # get variables
          f2 <- as.formula(paste0(wcols$response, " ~ -1 + ", wcols$covars))
          # calibration
          pred_vals <- model.matrix(f2, test_data) %*% coef(ridge_fit)[-1,] # multiple fits?

          # here we do not do it for all the metrics, should we? ??????????????
          d <- apply(pred_vals, 2, function(x = x, y = y, strat = strat){
            metric_fun(data.frame(x = x, y = y, strat = strat), errors=F)},
            y = test_data[[wcols$response]], strat = rep(1, nrow(test_data)))

          # coefficients
          if(metric == "coxnet.deviance") opt_fun <- which.min else opt_fun <- which.max
          coef_ridge <- matrix(coef(ridge_fit, s = ridge_fit$lambda[opt_fun(d)])[-1, , drop=FALSE]) # coefficients
          rownames(coef_ridge) <- rownames(coef(ridge_fit))[-1]

          #---- prepation to standardize coefs
          if(standardize == "internal" & tolower(method[1]) != "ridge") {

            if(verbose) print("Standardizing coefs...")

            # covariates summary
            all_vars <- all.vars(f2)[-1]
            classes <- sapply(data[,all_vars], class)
            # numeric variables
            data_summary_num <- as.data.frame(apply(na.omit(as.matrix(data[,all_vars[classes == "numeric"]])), 2, data_summary))
            # character variables - use mode
            data_summary_ch <- as.data.frame(apply(na.omit(as.matrix(data[,all_vars[classes != "numeric"]])), 2, data_summary_char))
            names(data_summary_ch) <- all_vars[classes != "numeric"]
            if(nrow(data_summary_ch) > 0) {
              dat_summ <- cbind(data_summary_num, data_summary_ch)[order(c(which(classes == "numeric"), which(classes != "numeric")))]
            } else {
              dat_summ <- data_summary_num
            }

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

        }

        # print(method[1])

        # if Adaptive Lasso
        if(tolower(method[1]) == "adaptivelasso") {

          if(verbose) print("Fitting AdaptiveLasso...")

          penalty.factor <- 1/(abs(coef_ridge)**gamma)
          penalty.factor[penalty.factor == Inf] <- 999999999 # If there is any infinite coefficient

        } else {

          # if OneZOI
          if(tolower(method[1]) == "onezoi-adaptivelasso" | tolower(method[1]) == "oz-adaptivelasso") {

            if(verbose) print("Fitting One-ZOI AdaptiveLasso...")

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
            for(jj in zoi_terms) {
              vals <- penalty.factor[mm_is_zoi == 1 & mm_predictor_vars == jj]
              # keep only the minimum
              vals[vals > min(vals, na.rm = TRUE)] <- Inf
              penalty.factor[mm_is_zoi == 1 & mm_predictor_vars == jj] <- vals
            }

            penalty.factor[penalty.factor == Inf] <- 999999999 # If there is any infinite coefficient
            # rownames(penalty.factor) <- mm_predictor_vars

          } else {

            # if grouped
            if(tolower(method[1]) == "grouped-adaptivelasso" | tolower(method[1]) == "g-adaptivelasso") {

              if(verbose) print("Fitting Grouped AdaptiveLasso...")

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
              for(jj in zoi_terms) {
                vals <- coef_ridge[mm_is_zoi == 1 & mm_predictor_vars == jj]
                # sum relative to the vatiation in the group
                vals <- grouped_func(vals, phi_group = factor_grouped_lasso)**gamma
                penalty.factor[mm_is_zoi == 1 & mm_predictor_vars == jj] <- vals
              }

              penalty.factor[penalty.factor == Inf] <- 999999999 # If there is any infinite coefficient

            } else {

              # if grouped
              if(tolower(method[1]) == "hypothesisdriven-adaptivelasso" | tolower(method[1]) == "hd-adaptivelasso") {

                if(verbose) print("Fitting HypothesisDriven AdaptiveLasso...")

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
                # expected_sign <- ifelse(mm_is_zoi == 1, expected_sign_hypothesis, 0)
                # penalty.factor <- hypothesis_func(coef_ridge, expected_negative, phi_hyp = factor_hypothesis)**gamma
                penalty.factor <- ifelse(mm_is_zoi == 1, function_hypothesis(-1*expected_sign_hypothesis*coef_ridge)**gamma,
                                         1/(abs(coef_ridge)**gamma))

                penalty.factor[penalty.factor == Inf] <- 999999999

              } else {

                if(tolower(method[1]) != "ridge") {
                  stop("Please choose a valid method for the function.")
                }

              }
            }
          }
        }
      } else {

        if(verbose) print("Fitting Lasso...")

      }
    }
  }

  # perform penalized regression with glmnet
  if(tolower(method[1]) == "ridge") {

    # get ridge fit
    fit <- ridge_fit

  } else {

    # fit
    fit <- net_logit(f, train_data,
                     alpha = alpha,
                     penalty.factor = penalty.factor,
                     type.measure = "deviance",
                     standardize = std,
                     na.action = na.action,
                     ...)

  }

  # get variables
  f2 <- as.formula(paste0(wcols$response, " ~ -1 + ", wcols$covars)) # should we remove the intercept?

  #----
  # Variable selection step

  # Calibration with conditionalBoyce index (recalled metric here)
  # does not work fro cv.glmnet
  # predict.glmnet(fit, test_data, type = "response")??
  pred_vals <- model.matrix(f2, test_data) %*% coef(fit)[-1,] # multiple fits?

  # variable names
  var_names <- rownames(coef(fit))[-1] # variable names

  # register used penalty.factor
  penalty_factor_modified <- penalty.factor
  if(!is.null(penalty_factor_modified)) names(penalty_factor_modified) <- var_names

  # prepare SDs for unstardization
  # if(standardize != FALSE) {
  if(standardize == "external") {

    # get variable name for each term
    # term_labels <- attr(terms(result$formula_no_strata), "term.labels")
    term_labels <- var_names
    variables <- unlist(lapply(term_labels, function(term) all.vars(parse(text = term)[[1]])))
    # check if there are categorical variables
    cat_vars <- names(numeric_covs[!numeric_covs])
    if(length(cat_vars) > 0) {
      cat_var_order <- list()
      for(i in seq_along(cat_vars)) {
        cat_var_order[[i]] <- indexes <- grep(cat_vars[i], variables)
        variables[indexes] <- cat_vars[i]
      }
    }
    # unique set of variable names
    variables_names <- unique(variables)

    sds <- covs_mean_sd$sd
    # if there are categorical variables
    if(length(cat_vars) > 0) {
      for(i in seq_along(cat_vars)) {
        sds <- append(sds, 1, after = cat_var_order[[i]][1]-1)
      }
    }

    # Create a named vector
    sd_vars <- setNames(sds, variables_names)  # Named vector: x -> 10, z -> 20

    # repeat values based on occurrences in the formula
    sds_all_terms <- sd_vars[variables]

  }

  # compute optimal score for multiple selected metrics
  metrics_evaluated <- list()
  for(mt in metrics_evaluate) {

    # get metric function
    mt_fun <- getFromNamespace(mt, ns = "oneimpact")

    # set min or max as optim function
    if(mt == "coxnet.deviance") opt_fun <- which.min else opt_fun <- which.max

    # compute the scoring metric using the test data
    d <- apply(pred_vals, 2, function(x = x, y = y, strat = strat){
      mt_fun(data.frame(x = x, y = y, strat = strat), errors=F)},
      y = test_data[[wcols$response]], strat = rep(1, nrow(test_data)))

    # save results
    metrics_evaluated[[mt]]$metric <- mt
    metrics_evaluated[[mt]]$test_scores <- d
    metrics_evaluated[[mt]]$opt_fun <- opt_fun
    metrics_evaluated[[mt]]$lambda_opt <- fit$lambda[opt_fun(d)]

    # print
    if(verbose) {
      print(metrics_evaluated[[mt]]$metric)
      print(metrics_evaluated[[mt]]$lambda_opt)
      plot(fit$lambda, metrics_evaluated[[mt]]$test_scores); abline(v = metrics_evaluated[[mt]]$lambda_opt)
      plot(fit, xvar = "lambda"); abline(v = metrics_evaluated[[mt]]$lambda_opt)
      pie(abs(coef(fit, s = metrics_evaluated[[mt]]$lambda_opt)[,1]), labels = rownames(fit$beta))
    }

    # coefs
    coef <- matrix(coef(fit, s = fit$lambda[opt_fun(d)])[-1,, drop=FALSE]) # coefficients
    rownames(coef) <- rownames(coef(fit, s = fit$lambda[opt_fun(d)]))[-1]

    # get predicted values based on the training, testing, and validation data
    train_pred_vals <- model.matrix(f2, train_data) %*% coef
    test_pred_vals <- model.matrix(f2, test_data) %*% coef
    val_pred_vals <- model.matrix(f2, validate_data) %*% coef

    # if the standardization was not done before calling the function
    if(standardize != FALSE) {

      if(standardize == "external") {
        coef_std <- coef
        coef <- apply(coef, MARGIN = 2, oneimpact::from_std, sd = sds_all_terms)
      } else {
        # coef_std <- apply(coef, MARGIN = 2, oneimpact::from_std, sd = sds_all_terms)
        coef_std <- NULL
      }

    } else {

      coef_std <- NULL

    }

    # register coefficients (standardized and unstandardized)
    metrics_evaluated[[mt]]$coef <- coef
    metrics_evaluated[[mt]]$coef_std <- coef_std

    # if the metric is coxnet.deviance, use it for the tuning parameter but compute
    # the validation score using the Cindex
    if(mt == "coxnet.deviance") mt_fun <- oneimpact::Cindex

    # save results train and test scores
    metrics_evaluated[[mt]]$train_score <- mt_fun(data.frame(x = train_pred_vals[,1],
                                                             y = train_data[[wcols$response]],
                                                             strat = rep(1, nrow(train_data))))
    metrics_evaluated[[mt]]$test_score <- unname(d[opt_fun(d)])

    #----
    # Validation step
    val <- data.frame(x = val_pred_vals[,1],
                      y = validate_data[[wcols$response]],
                      strat = rep(1, nrow(validate_data)))

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
      metrics_evaluated[[mt]]$validation_score <- unlist(lapply(val2, mt_fun))

    } else {
      metrics_evaluated[[mt]]$validation_score <- mt_fun(val)
    }

    # Validation habitat only
    # no habitat validation score
    metrics_evaluated[[mt]]$habitat_validation_score <- NULL
  }

  # initiate results object
  results <- list()

  # register parameters and all resuts
  results$parms <- parms
  results$glmnet_fit <- fit

  # all metrics evaluated
  results$metrics_evaluated <- metrics_evaluated
  results$var_names <- var_names # variable names
  results$penalty_factor_modified <- penalty_factor_modified
  # Add info about the covariates - type
  results$numeric_covs <- numeric_covs

  # standarized means and sd
  if(standardize != FALSE) {
    results$covariate_mean_sd <- covs_mean_sd
  } else {
    results$covariate_mean_sd <- NULL
  }

  # lambda and coefs for the selected metric
  results$metric <- metric
  results$alpha <- alpha
  results$lambda <- metrics_evaluated[[metric]]$lambda_opt

  results$coef <- metrics_evaluated[[metric]]$coef
  results$coef_std <- metrics_evaluated[[metric]]$coef_std

  #in case there are missing variables, add missing coefs
  if(any_missing_vars) {

    if(replace_missing_NA) {

      # get variables form the original formula
      f3 <- as.formula(paste0(wcols$response, " ~ -1 + ",
                              extract_response_strata(f_original, covars = TRUE)$covars)) # should we remove the intercept?
      # create original model matrix
      mm <- colnames(model.matrix(f3, train_data))

      # re-write coefs with NAs
      coef <- matrix(nrow = length(mm), ncol = 1)
      rownames(coef) <- mm
      # fill in the estimated coefficients
      coef[match(rownames(results$coef), mm)] <- results$coef
      results$coef <- coef # replace

      # standardized coefficients
      if(!is.null(coef_std)) {
        # re-write coefs_std with NAs
        coef_std <- matrix(nrow = length(mm), ncol = 1)
        rownames(coef_std) <- mm
        # fill in the
        coef_std[match(rownames(results$coef_std), mm)] <- results$coef_std
        results$coef_std <- coef_std
      }

    }
  }

  results$train_score <- metrics_evaluated[[metric]]$train_score
  results$test_score <- metrics_evaluated[[metric]]$test_score
  results$validation_score <- metrics_evaluated[[metric]]$validation_score
  results$validation_score_avg <- mean(results$validation_score, na.rm = TRUE)
  results$habitat_validation_score <- NULL
  results$habitat_validation_score_avg <- NULL

  # lambda and coefs for all metrics
  results$lambdas <- sapply(metrics_evaluated, function(x) x$lambda_opt)
  results$coefs_all <- sapply(metrics_evaluated, function(x) x$coef)
  results$coefs_std_all <- sapply(metrics_evaluated, function(x) x$coef_std)
  rownames(results$coefs_all) <- var_names
  if(standardize == "external") rownames(results$coefs_std_all) <- var_names
  results$train_score_all <- sapply(metrics_evaluated, function(x) x$train_score)
  results$test_score_all <- sapply(metrics_evaluated, function(x) x$test_score)
  results$validation_score_all <- sapply(metrics_evaluated, function(x) x$validation_score)
  results$habitat_validation_score_all <- NULL

  # whether to save the results externally
  if (!is.null(out_dir_file)){
    names_out <- oneimpact:::pretty_seq(1:999)[i]
    saveRDS(results, file = paste0(out_dir_file, "_", names_out, ".rds"))
  } else {
    return(results)
  }
}

#' @rdname fit_net_functions
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
#' @name bag_fit_net_functions
#' @export
bag_fit_net_logit <- function(f, data,
                              samples,
                              metric = c(AUC, conditionalBoyce, conditionalSomersD, conditionalAUC)[[1]],
                              method = c("Lasso", "Ridge", "AdaptiveLasso", "DistanceDecayLasso", "ElasticNet")[1],
                              standardize = c("internal", "external", FALSE)[1],
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
  numeric_covs <- (sapply(data_covs, is.numeric) == TRUE)
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
    ### warning if the is any cov with sd = 0, remove it or bring an error
    # if(covs_mean_sd)
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
  results$samples <- samples
  results$standardize <- standardize
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
                        predictor_table = predictor_table,
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
                                            predictor_table = predictor_table,
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

#' Function to set penalties according to hypotheses
#'
#' @param coefs Vector of estimated coefficients through a Ridge regression.
#' @param expectation Expected/Hypothetical signal of the coefficient for a
#' given covariate. Should take value `hyp = -1` or `hyp = +1` depending
#' on the expected signal of the coefficient. If `expectation = 0`, no change is
#' made and the penalty if `penalty.factor = 1/abs(coef)**gamma`.
#' @param phi_hyp Additional penalty constant for the hypothesis-based penalties.
#' A value in the interval [1, Inf] where 1 is no additional penalty and
#' higher values correspond to higher penalties when
#'
#' @examples
#' # set coefficients
#' coefs <- c(-1, -0.5, -0.1, 0.8, 0.3, -0.1)
#' expected_sign <- -1
#' hypothesis_func(coefs)
#'
#' x <- seq(-2, 2, length.out = 101)
#' plot(x, exp(x), ylab = "Penalty", xlab = "Coefficient")
#' plot(x, hypothesis_func(x), ylab = "Penalty", xlab = "Coefficient")
#' plot(x, hypothesis_func(x, phi_hyp = 50), ylab = "Penalty", xlab = "Coefficient")
#' plot(x, exp(x)*hypothesis_func(x, phi_hyp = 10), ylab = "Penalty", xlab = "Coefficient")
#'
#' @keywords internal
#' @export
hypothesis_func <- function(coefs, expectation = -1, phi_hyp = 1) {

  1/ifelse(coefs/expectation == abs(coefs) | expectation == 0, abs(coefs), 1/phi_hyp * abs(coefs))
  #ifelse(coefs/hyp == abs(coefs), 1/abs(coefs), 1/abs(coefs) - 1/abs(coefs) * phi_hyp)

}

#' @param phi_group Additional penalty constant for the group-based penalties.
#' A value in the interval [0, Inf] where 0 is no additional penalty and
#' higher values correspond to higher penalties.
#'
#' @rdname fit_net_functions
#' @export
grouped_func <- function(coefs, phi_group = 0) {

  # coefs_sorted <- sort(coefs, decreasing = TRUE)
  max_coef <- max(abs(coefs))
  min_coef <- min(abs(coefs))
  1/abs(coefs) + phi_group * (max_coef - abs(coefs))/(max_coef - min_coef)
  #ifelse(coefs/hyp == abs(coefs), 1/abs(coefs), 1/abs(coefs) - 1/abs(coefs) * phi_hyp)

}
