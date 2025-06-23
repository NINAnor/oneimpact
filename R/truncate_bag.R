#' Truncate bag to avoid weirdness in the model
#'
#' This function identifies sources of weirdness in the models in a bag and
#' uses them to remove variables the produce such weirdness. Sources of
#' weirdness might be coefficients with signs opposite to one's hypothesis,
#' response curves crossing zero, or response curves with multiple inflection
#' points, for instance.
#'
#' Currently, this function is applied only to terms corresponding to zones
#' of influence (ZOI). Importantly, this function does not re-fit the
#' model, but only sets coefficients to zero for all ZOI terms above the
#' radius in which a certain weirdness is identified.
#'
#' @param x Bag.
#' @param data `[data.frame]` \cr The original, complete data used for model fitting.
#' @param measure `[character(1)="cross"]{"coef_sign", "cross"}` \cr Measure used
#' to quantify "weirdness" in the model or coefficients, based on the coefficients and
#' the response plots for each type of covariate with zone of influence in a model.
#' It can be one of these:
#' - `"coef_sign"`: The measure is based on the minimum ZOI radius for which the
#' sign is opposite to the ecologically expected sign;
#' - `"cross"`: default. The measure is based on the minimum distance at which a
#' reponse curve crosses zero.
#' @param criterion `[character(1)="first_coef"]{"min", "first_coef"}` \cr Criterion
#' used to truncate the curves/coefficients, for each type of ZOI variable.
#' Only applicable for `measure = "cross"`.
#' If `criterion = "first_coef"` (default), the coefficients are set to zero
#' starting from the first coefficient whose sign is opposite to the expected sign,
#' which is larger than the distance at which the response plot crosses zero
#' for that ZOI variable.
#' If `criterion = "min"`, the coefficients are set to zero for all terms whose
#' radius is larger than the distance at which the response plot crosses zero
#' for that ZOI variable (regardless of the coefficient signs).
#' @param wmean `[logical(1)=TRUE]` \cr Whether the truncation should be based on the weighted mean
#' coefficients and response plots (default, if `wmean = TRUE`) or on each individual
#' model coefficient and response plots (if `wmean = FALSE`).
#' @param expected_sign `[numeric(1)=-1]` \cr Expected sign of the coefficient. Either -1 (negative),
#' +1 (positive), or 0 (no effect).
#' @param reassess `[logical(1)=TRUE]` \cr Should the model be reassessed after
#' truncation, with fit, calibration, and validation scores re-computed?
#' Default is `TRUE`.
#' @param ... \cr Other parameters used in [oneimpact::weirdness()].
#'
#' @example examples/truncate_bag_example.R
#'
#' @seealso [oneimpact::weirdness()]
#'
#' @export
truncate_bag <- function(x,
                         data,
                         measure = c("coef_sign", "cross")[2],
                         criterion = c("min", "first_coef")[2],
                         wmean = TRUE,
                         expected_sign = -1,
                         reassess = TRUE,
                         ...) {

  # compute weirdness
  weird <- weirdness(x = x,
                     data = data,
                     wmean = wmean,
                     expected_sign = expected_sign,
                     ...)

  # modified bag
  new_bag <- x

  # get ZOI variables and terms from predictor table
  pred_table <- x$parms$predictor_table

  # modifiy coefficients
  if(measure == "coef_sign") {
    weird$coef_sign_radii
    stop("To be implemented.")
  } else {

    if(measure == "cross") {

      truncate_points <- weird$where_crosses
      for(i in seq_along(truncate_points)) {
        var <- truncate_points[[i]]
        var_name <- names(truncate_points)[i]

        if(wmean) {
          if(length(var) > 0) {
            # get terms with radius higher than the crossing
            terms_to_zero <- pred_table |>
              dplyr::filter(grepl(var_name, term_zoi),
                            zoi_radius >= min(var)) |>
              dplyr::pull(term_zoi)
            if(criterion == "first_coef") {
              coefs_term <- new_bag$coef[rownames(new_bag$coef) %in% terms_to_zero,] %*% new_bag$weights
              first_coefs_term_againt_expected <- min(weirdness(coefs_term, expected_sign = expected_sign, which_coef_sign = "index"))
              if(length(first_coefs_term_againt_expected) > 0) {
                terms_to_zero <- terms_to_zero[first_coefs_term_againt_expected:length(terms_to_zero)]
              }
            }
            new_bag$coef[rownames(new_bag$coef) %in% terms_to_zero,] <- 0

          }
        } else {

          col_weight_positive <- which(x$weights > 0)

          for(vv in seq_along(var)) {
            var_vv <- var[[vv]]
            var_vv_name <- names(var)[vv]
            if(length(var_vv) > 0) {
              # get terms with radius higher than the crossing
              terms_to_zero <- pred_table |>
                dplyr::filter(grepl(var_name, term_zoi),
                              zoi_radius >= min(var_vv)) |>
                dplyr::pull(term_zoi)
              if(criterion == "first_coef") {
                coefs_term <- new_bag$coef[rownames(new_bag$coef) %in% terms_to_zero,col_weight_positive[vv]]
                first_coefs_term_againt_expected <- min(weirdness(coefs_term, expected_sign = expected_sign, which_coef_sign = "index"))
                if(length(first_coefs_term_againt_expected) > 0) {
                  terms_to_zero <- terms_to_zero[first_coefs_term_againt_expected:length(terms_to_zero)]
                }
              }
              new_bag$coef[rownames(new_bag$coef) %in% terms_to_zero, col_weight_positive[vv]] <- 0
            }
          }
        }

      }
    }
  }

  if(reassess) {

    # get variables
    f <- new_bag$formula

    # get variables
    wcols <- extract_response_strata(f2, covars = TRUE)

    # case
    case <- wcols$response
    # get variable referring to strata
    strat <- wcols$strata

    # check if the model is logit or clogit
    call <- if(strat == "") "fit_net_logit" else "fit_net_clogit"

    # filter out NAs in the response variable
    if(anyNA(data[[wcols$response]])) {
      data <- data[!is.na(data[[case]]),]
    }

    # relevant columns
    all_vars <- all.vars(f)

    # in any case, we keep their mean and sd which might be useful
    if(call == "fit_net_logit") {
      all_covars <- all_vars[-1]
    } else {
      all_covars <- grep(wcols$strata, all_vars[-1], invert = TRUE, value = TRUE)
    }

    # get predictors
    data_covs <- data[, all_covars]
    # select numeric predictors to be standardized
    numeric_covs <- new_bag$numeric_covs

    #----
    # Evaluate predictions for each

    # register/use standardized covariates
    if(new_bag$standardize == "external") {
      # standadize covs
      # get standardization parms
      data_covs_num <- data_covs[, numeric_covs]
      # standardize
      data_covs_num_std <- lapply(1:ncol(data_covs_num), function(i) scale(data_covs_num[,i]))

      # merge standardized predictors with non numeric predictors
      data_covs_std <- cbind(data_covs[, !numeric_covs], data.frame(do.call("cbind", data_covs_num_std)))
      data_covs_std <- data_covs_std[,order(c(which(!numeric_covs), which(numeric_covs)))]
      colnames(data_covs_std) <- colnames(data_covs)
      if(call == "fit_net_clogit") {
        data <- cbind(data[wcols$response], data[wcols$strata], data_covs_std)
      } else {
        data <- cbind(data[wcols$response], data_covs_std)
      }

    } else {
      # or just use the original data
      data <- data[, all_vars]
    }

    # samples
    samples <- new_bag$samples

    #----------------------
    # loop over samples
    n_samp <- new_bag$n_no_errors

    i <- 1
    for(i in seq_len(n_samp)) {

      if(call == "fit_net_clogit") {
        # separate data for fitting, calibration, and validation
        if(is.null(new_bag$samples$sp_strat_id)) {
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
        }
        if(anyNA(test_data)) {
          n_bef <- nrow(test_data)
          test_data <- filter_na_strata(f, na.omit(test_data))
          nNA <- n_bef - nrow(test_data)
        }
        if(anyNA(validate_data)) {
          n_bef <- nrow(validate_data)
          validate_data <- filter_na_strata(f, na.omit(validate_data))
          nNA <- n_bef - nrow(validate_data)
        }


      } else {
        # separate data for fitting, calibration, and validation
        train_data  <- data[samples$train[[i]], all_vars]
        test_data <- data[samples$test[[i]], all_vars]
        validate_data <- data[samples$validate[[i]], all_vars]

        # check NAs
        if(anyNA(train_data)) {
          train_data <- na.omit(train_data)
          nNA <- length(na.action(train_data))
        }
        if(anyNA(test_data)) {
          test_data <- na.omit(test_data)
          nNA <- length(na.action(test_data))
        }
        if(anyNA(validate_data)) {
          validate_data <- na.omit(validate_data)
          nNA <- length(na.action(validate_data))
        }

      }

      # get variables
      f2 <- new_bag$formula_no_strata

      # only single metric
      metrics_evaluate <- new_bag$metrics
      # compute optimal score for multiple selected metrics
      metrics_evaluated <- list()
      mt <- metrics_evaluate[1]
      for(mt in metrics_evaluate) {

        # get metric function
        mt_fun <- getFromNamespace(mt, ns = "oneimpact")

        # set min or max as optim function
        if(mt == "coxnet.deviance") opt_fun <- which.min else opt_fun <- which.max

        # coefs
        coef <- new_bag$coef
        coef_std <- new_bag$coef_std

        # get predicted values based on the training, testing, and validation data
        train_pred_vals <- model.matrix(f2, train_data) %*% coef
        test_pred_vals <- model.matrix(f2, test_data) %*% coef
        val_pred_vals <- model.matrix(f2, validate_data) %*% coef

        if(call == "fit_net_clogit") {
          if(new_bag$parms$kernel_vars[1] != "") {
            pred_vals_kernel <- kernel_prediction(f, validate_data,
                                                  kernel_vars = kernel_vars,
                                                  coefs = coef[,1])
          }
        }

        # if the metric is coxnet.deviance, use it for the tuning parameter but compute
        # the validation score using the Cindex
        if(mt == "coxnet.deviance") mt_fun <- oneimpact::Cindex

        if(call == "fit_net_logit") {
          wcols$strata <- "strat"
          train_data[[wcols$strata]] <- 1
          test_data[[wcols$strata]] <- 1
          validate_data[[wcols$strata]] <- 1
        }

        # save new_bags train and test scores
        train_score <- mt_fun(data.frame(x = train_pred_vals[,1],
                                         y = train_data[[wcols$response]],
                                         strat = train_data[[wcols$strata]]))
        test_score <- mt_fun(data.frame(x = test_pred_vals[,1],
                                        y = test_data[[wcols$response]],
                                        strat = test_data[[wcols$strata]]))

        #----
        # Validation step
        val <- data.frame(x = val_pred_vals[,1],
                          y = validate_data[[wcols$response]],
                          strat = validate_data[[wcols$strata]])

        if(call == "fit_net_clogit") {

          val <- merge(val, data.frame(strat = samples$sp_strat_id, blockH0=samples$blockH0),
                       by = "strat", all.x = T, all.y = F)

          if(!is.null(samples$blockH0)) {

            val2 <- split(val, val$blockH0)

            if(length(val2) == 0) {
              if(is.null(samples$sp_strat_id)) {
                val2 <- split(val, samples$blockH0[samples$validate[[i]]])
              } else {
                val2 <- split(val, samples$validate[[i]])
              }
            }
            validation_score <- unlist(lapply(val2, mt_fun))

          } else {

            validation_score <- mt_fun(val)
          }

          # Validation habitat only
          if(new_bag$parms$kernel_vars[1] != "") {

            pred_vals_habitat <- val_pred_vals - pred_vals_kernel # does it make sense??
            hab <- data.frame(x = pred_vals_habitat[,1],
                              y = validate_data[[wcols$response]],
                              strat = validate_data[[wcols$strata]])
            hab <- merge(hab, data.frame(strat = samples$sp_strat_id, blockH0=samples$blockH0),
                         by = "strat", all.x = T, all.y = F)

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
              habitat_validation_score <- unlist(lapply(hab2, mt_fun))

            } else {
              habitat_validation_score <- mt_fun(hab)
            }

          } else {

            # if there is not movement kernel terms, NULL
            habitat_validation_score <- NULL

          }

        } else {
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
            validation_score <- unlist(lapply(val2, mt_fun))

          } else {
            validation_score <- mt_fun(val)
          }

          # Validation habitat only
          # no habitat validation score
          habitat_validation_score <- NULL
        }

      }

      new_bag$fit_score[i] <- train_score
      new_bag$calibration_score[i] <- test_score
      new_bag$validation_score[i] <- validation_score
      if(!is.null(new_bag$habitat_validation_score)) new_bag$habitat_validation_score[i] <- habitat_validation_score

    }

    ###
    # Do we need to reassess weights???

    weights_matrix <- matrix(rep(new_bag$weights, each = nrow(coef)), nrow = nrow(coef))
    new_bag$wcoef <- coef * as.vector(weights_matrix) # each model
    new_bag$wcoef_std <- coef_std * as.vector(weights_matrix) # each model

    # summary validation scores
    if(nrow(new_bag$validation_score) > 1) {
      new_bag$validation_score_summary <- data.frame(min = apply(new_bag$validation_score, 2, min,na.rm = TRUE),
                                                    median = apply(new_bag$validation_score, 2, median,na.rm = TRUE),
                                                    mean = apply(new_bag$validation_score, 2, mean,na.rm = TRUE),
                                                    max = apply(new_bag$validation_score, 2, max,na.rm = TRUE))
    } else {
      new_bag$validation_score_summary <- data.frame(median = apply(new_bag$validation_score,2,mean,na.rm = TRUE))
    }

    # weighted validation
    new_bag$weighted_validation_score <-  new_bag$validation_score %*% new_bag$weights
    colnames(new_bag$weighted_validation_score) <- "weighted_validation_score"

    if(nrow(new_bag$validation_score) > 1) {

      new_bag$weighted_validation_score_summary <- cbind(
        apply(new_bag$validation_score, 2, min, na.rm = TRUE) %*% new_bag$weights,
        apply(new_bag$validation_score, 2, median,na.rm = TRUE) %*% new_bag$weights,
        apply(new_bag$validation_score, 2, mean,na.rm = TRUE) %*% new_bag$weights,
        apply(new_bag$validation_score, 2, max,na.rm = TRUE) %*% new_bag$weights)
      colnames(new_bag$weighted_validation_score_summary) <- c("min", "median", "mean", "max")
    } else {
      new_bag$weighted_validation_score_summary <-
        apply(new_bag$validation_score,2,min,na.rm = TRUE) %*% new_bag$weights
      colnames(new_bag$weighted_validation_score_summary) <- "mean"
    }

  } # end of reassessment

  new_bag
}

