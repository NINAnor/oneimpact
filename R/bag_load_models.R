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
bag_load_models <- function(data,
                            load_models_path = ".",
                            load_models_pattern = NULL,
                            names_from_file = FALSE,
                            name_from_file_pattern = "Resample",
                            verbose = FALSE) {

  # initiate results object
  results <- list()

  #---
  # if the models were already run, read them
  model_files <- list.files(path = load_models_path, pattern = load_models_pattern,
                            full.names = TRUE) |>
    grep(pattern = ".rds$", value = TRUE)

  fitted_list <- list()
  for(i in seq_along(model_files)) {
    if(verbose) print(paste0("Loading model ", i, "/", length(model_files), "..."))
    fitted_list[[i]] <- readRDS(model_files[i])
  }

  #---
  # get number of models
  results$n <- length(fitted_list[[1]]$parms$samples$train)

  # check number of files
  if(length(model_files) != results$n)
    warning(paste0("Warning: there should be ", results$n, " models, but we found ", length(model_files), " files. Please check."))

  # add errors to the others - flag
  if(names_from_file) {
    model_number <- strsplit(model_files, split = name_from_file_pattern) |>
      sapply(function(x) x[2]) |>
      sapply(function(x) as.numeric(gsub("\\D", "", x)))
    names_out <- oneimpact::pretty_seq(1:999)[unname(model_number)]
    names(fitted_list) <- names_out
  } else {
    names(fitted_list) <- names(fitted_list[[1]]$parms$samples$train)
  }

  # get other parms
  results$formula <- fitted_list[[1]]$parms$f
  results$method <- fitted_list[[1]]$parms$method
  results$metric <- fitted_list[[1]]$parms$metric
  results$samples <- fitted_list[[1]]$parms$samples
  results$standardize <- fitted_list[[1]]$parms$standardize

  # set list of models
  results$models <- fitted_list

  results
}
