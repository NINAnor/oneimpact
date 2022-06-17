#' Evaluation of single covariate effects at multiple scales for multiple covariates
#'
#' @param covariates `[character]` \cr List of names of the covariates to be tested,
#' or a pattern within their names that differentiate among them.
#'
#' @export
multifit_single_multivar <- function(mod,
                                     covariates,
                                     data,
                                     formula = NULL,
                                     args = NULL,
                                     criterion = "AIC",
                                     corr_check = TRUE,
                                     corr_criterion = c("usdm::vifcor",
                                                        "usdm::vifstep")[1],
                                     corr_threshold = 0.7,
                                     print_best = 5, ...) {


  multifits <- list()
  model_comparison <- list()
  best_model <- list()
  model_elements <- list()

  # For each covariate
  for(i in 1:length(covariates)) {

    # covariate columns in the input dataset
    cov_names <- colnames(data) %>%
      grep(pattern = covariates[[i]], value = T)

    # run multifit_single
    multifit_var <- multifit_single(mod = mod,
                                    multief = cov_names,
                                    data = data,
                                    formula = formula,
                                    args = args,
                                    criterion = criterion,
                                    ...)

    multifits[[i]] <- multifit_var

    # summarize results
    model_comparison_df <- multifit_var$summary %>%
      tibble::rowid_to_column()

    # sort
    if(criterion %in% c("AIC", "BIC")) {
      if(criterion == "AIC") {
        model_comparison_df <- model_comparison_df %>%
          # dplyr::arrange(!!dplyr::sym(criterion)) %>%
          dplyr::arrange(AIC) %>%
          dplyr::mutate(dAIC = AIC - AIC[1])
      }
      if(criterion == "BIC") {
        model_comparison_df <- model_comparison_df %>%
          # dplyr::arrange(!!dplyr::sym(criterion)) %>%
          dplyr::arrange(BIC) %>%
        dplyr::mutate(dBIC = BIC - BIC[1])
      }
    } else {
      model_comparison_df <- model_comparison_df %>%
        dplyr::arrange(dplyr::desc(!!dplyr::sym(criterion)))
    }

    # get coefficients
    coefs <- purrr::map2(
      multifit_var$models,
      multifit_var$lands_summary$spatial_scale,
      function(x, y) {
        coefs <- x$coefficients
        unname(coefs[match(y, names(coefs))])
      }) %>% unlist()

    # get correlation between variables
    if(corr_check) {
      corr <- purrr::map2(
        multifit_var$models,
        multifit_var$summary$AIC,
        function(x, y) {
          if(!is.na(y)) {
            m_vars <- all.vars(x$formula)[-1]
            m_vars <- names(dplyr::select(data[, m_vars], where(is.numeric)))
            corr_m <- do.call(eval(parse(text = corr_criterion)), list(x = data[, m_vars], th = corr_threshold))
          } else {
            corr_m <- NA
          }
          corr_m
        })
    } else
      corr <- NULL

    # get model formula
    formulas <- lapply(multifit_var$models, function(x) x$formula)

    # add formula, coefficients, correlation, warnings, and messages
    model_comparison[[i]] <- model_comparison_df %>%
      tibble::as_tibble() %>%
      dplyr::mutate(beta = coefs[model_comparison_df$rowid])

    model_elements[[i]] <- list(corr = corr[model_comparison_df$rowid],
                                formula = formulas[model_comparison_df$rowid],
                                warnings = multifit_var$warnings[model_comparison_df$rowid],
                                messages = multifit_var$messages[model_comparison_df$rowid])

    # get best model
    best_model[[i]] <- multifit_var$models[[model_comparison_df$rowid[1]]]
  }

  # top models of each variable in a single
  model_comparison_all <- lapply(seq_along(model_comparison), function(i) {
    model_comparison[[i]][1:print_best,] %>%
      dplyr::mutate(covariate = covariates[i],
                    rank = 1:print_best) %>%
      dplyr::relocate(rank) %>%
      dplyr::relocate(covariate)
    }) %>%
    dplyr::bind_rows() %>%
    tibble::as_tibble()

  # rename elements
  names(multifits) <- covariates
  names(model_comparison) <- covariates
  names(best_model) <- covariates
  names(model_elements) <- covariates

  # return all lists
  list(fits = multifits, model_comparison = model_comparison,
       best_model = best_model,
       model_elements = model_elements,
       model_comparison_all = model_comparison_all)
}


grep
