#' Evaluation of multiple covariate effects at multiple scales
#'
#' @export
multifit_multiscale <- function(mod,
                                mod_comparison_df,
                                data,
                                formula = NULL,
                                weights = NULL,
                                args = NULL,
                                max_models = 10000,
                                corr_check = TRUE,
                                corr_criterion = c("usdm::vifcor", "usdm::vifstep")[1],
                                corr_threshold = 0.7,
                                progress = FALSE,
                                verbose = FALSE,
                                ...) {

  # take variables and their combination
  covs <- unique(mod_comparison_df$covariate)

  # check number of multivariate models to be run
  n_per_cov <- purrr::map_dbl(covs, function(i)
    nrow(mod_comparison_df[mod_comparison_df$covariate == i,]))
  n_models <- prod(n_per_cov)
  if(n_models > max_models)
    stop(paste0("The number of models to be tested is higher than ", max_models, " (n = ", n_models,"). Please increase 'max_models' or rethink your approach."))

  # list of covariates
  cov_list <- purrr::map(covs, function(i)
    mod_comparison_df %>%
      dplyr::filter(covariate == i) %>%
      pull(multief) %>%
      as.character())
  names(cov_list) <- covs

  # list of parms
  parms <- expand.grid(cov_list) %>%
    tibble::as_tibble() %>%
    dplyr::mutate_all(as.character)

  # output result table
  out_df <- parms %>%
    dplyr::mutate(rank = as.numeric(NA),
                  AIC = as.numeric(NA),
                  dAIC = as.numeric(NA))

  # list of models
  models <- list()
  excluded_vars <- list()
  warns <- list()
  formulas <- list()

  # check args
  if(!is.null(args))
    initial.args <- paste(",", paste(args, collapse = ","), sep = "") else
    initial.args <- NULL

  if(progress) pb <- txtProgressBar(min = 0, max = nrow(out_df), style = 3)
  for (i in 1:nrow(out_df)) {

    # get variables
    vars <- unname(unlist(out_df[i,1:length(covs)]))

    # if correlation between variables should be tested
    # SO FAR: only checking among the infrastructure variables, not the other ones...
    if(corr_check) {
      # check correlated variables
      corr_m <- do.call(eval(parse(text = corr_criterion)), list(x = data[, vars], th = corr_threshold))
      # record variables to exclude
      exclude <- corr_m@excluded

      # update - keep only uncorrelated variables
      vars <- vars[!(vars %in% exclude)]
    } else {
      exclude <- NA
    }

    # update formula
    # add transformation to the parameters!!! e.g. scale
    ff <- as.character(formula)[c(2,1,3)]
    for(f in length(covs):1)
      ff[3] <- paste0(vars[[f]][1], " + ", ff[3])
      # update(ff, ~ . + scale(eval(parse(vars[[f]]))))
    # ff <- formula(paste(ff, collapse = " "))

    # run model
    running <- function(expr) {
      warns <- mess <- NULL
      warn_handler <- function(w) {
        warns <<- c(warns, list(w))
        invokeRestart("muffleWarning")
      }
      mess_handler <- function(m) {
        mess <<- c(mess, list(m))
        NULL
      }
      val <- suppressMessages(tryCatch(withCallingHandlers(expr, warning = warn_handler, message = mess_handler), error = function(e) e))
      out <- list(value = val, warnings = warns, messages = mess)
      return(out)
    }

    if(!is.null(weights))
      expression <- paste0(mod, "(", paste(ff, collapse = " "), ", data = ", deparse(substitute(data)),
                           ", weights = ", deparse(substitute(weights)), initial.args, ")") else
        expression <- paste0(mod, "(", paste(ff, collapse = " "), ", data = ", deparse(substitute(data)), initial.args, ")")

    if(verbose) print(expression)
    model <- try(running(eval(parse(text = expression))))

    # get models, formulas, warnings, excluded variables
    models[[i]] <- model$value
    excluded_vars[[i]] <- exclude
    formulas[[i]] <- try(model$value$formula)
    warns[[i]] <- model$warnings
    # get fit details
    out_df$AIC[i] <- try(AIC(model$value))

    if(progress) setTxtProgressBar(pb, i)
  }

  # errors in fitting
  which_error <- grep("Error", out_df$AIC)

  if(length(which_error) > 0) {
    out_df <- out_df %>%
      dplyr::slice(-which_error)
  }

  # organize results
  out_df <- out_df %>%
    tibble::rowid_to_column(var = "modelID") %>%
    dplyr::arrange(AIC) %>%
    dplyr::mutate(rank = 1:nrow(.),
                  dAIC = AIC - AIC[1],
                  relLL = exp(-0.5 * dAIC),
                  wAIC = relLL/sum(relLL)) %>%
    dplyr::select(-relLL)

  # return list with results and model fits
  list(model_comparison_df = out_df, best_model = models[[out_df[1,1][[1]]]], models = models, formulas = formulas, excluded_vars = excluded_vars, warnings = warns)
}
