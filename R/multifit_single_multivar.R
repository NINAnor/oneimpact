#' Evaluation of single covariate effects at multiple scales for multiple covariates
#'
#' @param covariates `[character]` \cr List of names of the covariates to be tested, or a pattern
#' within their names that differentiate among them.
#'
#' @export
multifit_single_multivar <- function(mod, covariates, data, formula = NULL, args = NULL, criterion = "AIC", print_best = 5, ...) {


  multifits <- list()
  model_comparison <- list()
  best_model <- list()

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

    # add model call, warnings, and messages
    calls <- lapply(multifit_var$models, function(x) x$call)

    model_comparison[[i]] <- model_comparison_df %>%
      tibble::as_tibble() %>%
      dplyr::mutate(call = calls[model_comparison_df$rowid],
                    warnings = multifit_var$warnings[model_comparison_df$rowid],
                    messages = multifit_var$messages[model_comparison_df$rowid])

    # get best model
    multifit_var$models[[model_comparison_df$rowid[1]]]
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

  best_model <- lapply()

  # rename elements
  names(multifits) <- covariates
  names(model_comparison) <- covariates

  # return all lists
  list(fits = multifits, model_comparison = model_comparison, best_model = best_model, model_comparison_all = model_comparison_all)
}
