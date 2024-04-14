#' Load a series of files output of fit_net_clogit models and put them on a bag
#'
#' @param f `[formula]` \cr Formula of the models fitted, with all possible candidate terms.
#' @param data `[data.frame,tibble]` \cr Complete data set analyzed.
#' @param samples `[list]` \cr List of samples with at least three elements: train, test,
#' and validate. Each elements might have several elements, each representing
#' the lines of `data` to be sampled for each resample. Typically, this is computed by
#' the function [oneimpact::create_resamples()].
#' @param load_models_path `[string="."]` \cr Path to the folder where the files
#' are saved.
#' @param load_models_pattern `[string="."]` \cr Pattern common to the file names,
#' to be used in `grep()` to find the files within the folder
#' `load_models_path`. It should be the pattern present in the argument
#' `out_dir_file` in [oneimpact::fit_net_clogit()] or [oneimpact::bag_fit_net_clogit()]
#' functions.
#' @param names_from_file `[logical(1)=TRUE]` \cr If `FALSE` (default), the names of
#' the resamples are taken from the `samples` parameter. If `TRUE`, they are taken
#' from the files names, instead. In this case, the string defined by the parameter
#' `name_from_file_pattern` is used to identify the number of each resample.
#' @param name_from_file_pattern `[string="Resample"]` \cr String used to separate the
#' file name and identify the number of the resample of each model in the bag, when
#' they are read from files.
#'
#' The parameters metric, standardize, and method should be same ones used to
#' fit the bag of models.
#'
#' @export
bag_load_models <- function(f, data, samples,
                            load_models_path = ".",
                            load_models_pattern = NULL,
                            names_from_file = FALSE,
                            name_from_file_pattern = "Resample",
                            metric = c(conditionalBoyce, conditionalSomersD, conditionalAUC)[[1]],
                            standardize = c("internal", "external", FALSE)[1],
                            method = c("Lasso", "Ridge", "AdaptiveLasso", "DecayAdaptiveLasso", "ElasticNet")[1],
                            verbose = FALSE) {

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

  # if the models were already run, read them
  model_files <- list.files(path = load_models_path, pattern = load_models_pattern,
                            full.names = TRUE) |>
    grep(pattern = ".rds$", value = TRUE)

  # initiate results object
  results <- list()
  results$n_samples <- length(samples$train)
  results$n <- length(model_files)
  results$formula <- f
  results$method <- method
  results$metric <- metric

  # standarized means and sd
  if(standardize == "external") {
    results$covariate_mean_sd <- covs_mean_sd
  } else {
    results$covariate_mean_sd <- NULL
  }

  # check number of files
  if(length(model_files) != results$n)
    warning(paste0("Warning: there should be ", results$n, " models, but we found ", length(model_files), " files. Please check."))

  fitted_list <- list()
  # for(i in 1:length(samples$train))
  for(i in 1:results$n) {
    if(verbose) print(paste0("Loading model ", i, "/", length(samples$train), "..."))
    fitted_list[[i]] <- readRDS(model_files[i])
  }

  # add errors to the others - flag
  if(names_from_file) {
    model_number <- strsplit(model_files, split = name_from_file_pattern) |>
      sapply(function(x) x[2]) |>
      sapply(function(x) as.numeric(gsub("\\D", "", x)))
    names_out <- oneimpact:::pretty_seq(1:999)[unname(model_number)]
    names(fitted_list) <- names_out
  } else {
    names(fitted_list) <- names(samples$train)
  }

  # define new class?
  results$models <- fitted_list

  ## TO DO
  # unstandardize coeffients if standarize = "external"

  # Add info about the covariates - type
  results$numeric_covs <- numeric_covs

  results
}
