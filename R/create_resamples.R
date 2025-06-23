#' Create samples for fitting, calibrating, and validating models
#'
#' The function creates (sub)samples of the data to be included for three different
#' model blocks: train (fitting), test (calibrate, tunning), and validate.
#' By default, samples are created through bootstrapping, i.e. with replacement. This means
#' data observations can be repeated within a given sample block, but observations included
#' in one block are necessarily excluded from the other blocks (e.g. observations selected
#' for validation will be absent from train and test blocks).
#''
#' Samples can be created at random (if `spat_strat = NULL`, default) or with spatial
#' stratification (spatial strata can be created with the function [oneimpact::spat_strat()].
#' In the latter case, train and test sets are spatially split, to allow for a more
#' thorough cross-validation to define the penalty parameter in the penalized regressions.
#' Also, samples might include a specific variable (with classes or groups) H0 to be used for
#' (block cross-)validation (if `colH0` is provided), but this is not a requirement.
#'
#' @param y `[vector]` \cr A vector of outcomes. It can be the response variable for
#' the data set of interest, or only the `case = 1` cases for conditional logistic (
#' step-selection) analyses.
#' @param times `[numeric(1)=10]` \cr The number of partitions or samples to be sampled.
#' @param p `[numeric(3)=c(0.4,0.2,0.2)]` \cr A 3 element numeric vector with the percentage of data that goes to
#' fitting/training (H1), testing (H2), and validation (H0). Values should be between 0 and 1 and should
#' not sum more than 1.
#' @param max_size_blockH0_validation `[numeric(1)=1000]` \cr Maximum size of the blocks H0 (e.g. population,
#' area, year) for validation block H0. Used to limit the number of observations in the validation set, to
#' avoid sampling too many observations of the block H0 levels with more observations,
#' for imbalanced data sets. To find out about meaningful values for this parameter,
#' use [oneimpact::explore_blocks_pre()] and [oneimpact::explore_blocks()].
#' @param max_size_blockH0_train `[numeric(1)=1000]` \cr Maximum size of the blocks H0 (e.g. population,
#' area, year) for training/fitting the model. Used to limit the number of observations in the train set,
#' to avoid sampling too many observations of the block H0 levels with more observations,
#' for imbalanced data sets. To find out about meaningful values for this parameter,
#' use [oneimpact::explore_blocks()].
#' @param max_size_blockH0_test `[numeric(1)=1000]` \cr Not implemented yet.
#' @param max_number_blocksH1_train `[numeric(1)=15]` \cr Maximum number of levels or blocks H1 to be used
#' for model fitting/training. This is only meaningful if there is spatial stratification
#' (i.e. if `sp_strat` in not `NULL`). To find out about meaningful values for this parameter,
#' use [oneimpact::explore_blocks()].
#' @param sp_strat `[data.frame]` \cr Default is `NULL`. If not `NULL`, the `data.frame` resulting
#' from [oneimpact::spat_strat()] should be provided here.
#' @param colH0 `[numeric,character,vector]` \cr Column number or name to define the IDs of the H0 level -
#' the one with ecological meaning, e.g. individual, population, or study area, used for validating the
#' predictions of the fitted model. If `sp_strat` is provided,
#' `colH0` is a string with the column name (or the column number) in the `sp_strat` table.
#' If `sp_strat = NULL`, `colH0` is a vector of H0 values with the same length as `y`.
#' If `colH0 = NULL` (Default), no H0 level is defined and there is no block cross-validation
#' in the bootstraped sets.
#' @param replace `[logical(1)=TRUE]` \cr Whether to perform the bootstrap sampling with or without
#' replacement (Default is `TRUE`).
#' @param H0setup Not implemented yet.
#'
#' @returns A list with lists for the sets for train, test, and validation, each of which
#' with the indices corresponding to the observations to be kept in each resample.
#' If `colkH0` is not `NULL`, a vector with the blockH0 which each observation pertains
#' to is also appended to the output.
#' If `spat_strat` is provided, a list of blocks H0 and possibly a list of strata
#' might also be provided.
#'
#' @examples
#' # random sampling, no validation block H0
#' y <- runif(200)
#' samples <- create_resamples(y, p = c(0.4, 0.2, 0.2), times = 5)
#' samples
#'
#' # with validation block H0
#' data(reindeer)
#' library(terra)
#' library(amt)
#'
#' # random sampling, with validation block H0
#' samples <- create_resamples(1:nrow(reindeer), times = 5,
#'                             p = c(0.2, 0.2, 0.2),
#'                             max_size_blockH0_validation = 1000,
#'                             colH0 = reindeer$original_animal_id)
#' samples
#'
#' # spatially stratified sampling, with validation block H0
#' spst <- spat_strat(reindeer, coords = c("x", "y"), colH0 = "original_animal_id",
#'                    all_cols = F)
#' samples <- create_resamples(1:nrow(reindeer), times = 5,
#'                             p = c(0.2, 0.2, 0.2),
#'                             max_number_blocksH1_train = 20,
#'                             sp_strat = spst,
#'                             colH0 = "blockH0")
#' samples
#' sum(is.na(samples$test[[1]]))
#' sapply(samples$train, function(x) sum(is.na(x)))
#' sapply(samples$test, function(x) sum(is.na(x)))
#'
#' # small number of blocks or too high p[1] might incur in errors
#' samples <- create_resamples(1:nrow(reindeer), times = 10,
#'                             max_number_blocksH1_train = 3,
#'                             sp_strat = spst,
#'                             colH0 = "blockH0")
#'
#' @export
create_resamples <- function (y, times = 10,
                              p = c(0.4, 0.2, 0.2),
                              max_size_blockH0_validation = 1000,
                              max_size_blockH0_train = 1000,
                              max_size_blockH0_test = 1000, # not implemented
                              max_number_blocksH1_train = 40,
                              sp_strat = NULL,
                              colH0 = NULL,
                              H0setup = c("LAO", "LOO")[1],
                              replace = TRUE) {
  if (inherits(y, "Surv"))
    y <- y[, "time"]

  # set number of points in each
  # here we can use different rules
  N <- floor(p * length(y))

  # set matrices
  # 3 matrices, for train, test, and validate
  indexes <- lapply(N, function(n) matrix(0, ncol = times, nrow = n))

  #-------------------------
  # random sampling
  if(is.null(sp_strat)) {

    print("Starting random sampling...")

    #-------------------------
    # random sampling, no validation level H0
    if(is.null(colH0)) {
      # sample
      # first validate, then train and test, mutually exclusive
      for(col in seq(times)) {

        # validate
        index_val <- seq(along = y)
        val <- sort(sample(index_val, size = nrow(indexes[[3]]), replace = replace))

        # train/fit
        index_train <- index_val |>
          setdiff(val)
        train <- sort(sample(index_train, size = nrow(indexes[[1]]), replace = replace))

        # test/calibrate
        index_test <- index_train |>
          setdiff(train)
        test <- sort(sample(index_test, size = nrow(indexes[[2]]), replace = replace))

        indexes[[1]][,col] <- train
        indexes[[2]][,col] <- test
        indexes[[3]][,col] <- val
      } # end of loop for times

    } else {

      #-------------------------
      # random sampling, with level H0

      # sample
      # first validate, then train and test, mutually exclusive
      for(col in seq(times)) {

        # validate
        index_val <- seq(along = y)

        # (sample validation block from populations (i.e. blockH0))
        # split by H0 level
        blocksH0 <- split(index_val, colH0)
        # LAO
        # sample pts in each pop/level - controlling for too imbalanced levels using max_size_blockH0_validation
        # size is also controlled by the proportion of the whole data set for validation, p[3]
        validation_set <- lapply(blocksH0, function(x){
          sample(x, size = min(max_size_blockH0_validation, length(x)), replace = replace)})
        # sample only the proportion set for validation
        val <- sort(sample(unname(unlist(validation_set)), size = N[3]))
        if(length(val) > N[3]) {
          val <- sort(sample(val, N[3]))
        } else {
          if(length(val) < N[3]) {
            warning(paste0("The size of the validation set is smaller than intended; you should check the values of p[3] (", p[3]
                           ,") and/or max_size_blockH0_validation (", max_size_blockH0_validation, ").",
                           "As a reference, the number of points assigned to validation is ", length(val), ", and the number of points wanted for validation is ", N[3], "."))
          }
        }

        # train/fit
        index_train <- index_val |>
          setdiff(val)
        train <- sort(sample(index_train, size = nrow(indexes[[1]]), replace = replace))

        # test/calibrate
        index_test <- index_train |>
          setdiff(train)
        test <- sort(sample(index_test, size = nrow(indexes[[2]]), replace = replace))

        indexes[[1]][,col] <- train
        indexes[[2]][,col] <- test
        indexes[[3]][,col] <- val
      } # end of loop for times
      # To implement, add condition for LOO and LAO
    } # end of condition for presence of H0
  } # end of condition for random sampling

  #-------------------------
  # spatial sampling, with level H0

  #-------------------------
  # spatially stratified sampling
  if(!is.null(sp_strat)) {

    print("Starting sampling with spatial stratification...")

    #-------------------------
    # spatially stratified sampling, no validation level H0
    if(is.null(colH0)) {
      stop("Spatial sampling without a given level H0 to be implemented.")
      # To implement, add condition for LOO and LAO
    } else {
      #-------------------------
      # spatially stratified sampling, with validation level H0

      # first validate, then train and test, mutually exclusive
      for(col in seq(times)) {

        #------------------
        # This is LAO, must set LOO in an alternative first block

        # validate
        # (sample validation block from populations (i.e. blockH0))
        # split by H0 level
        blocksH0 <- split(sp_strat, sp_strat[[colH0]])
        # LAO
        # sample pts in each pop/level - controlling for too imbalanced levels using max_size_blockH0_validation
        # size is also controlled by the proportion of the whole data set for validation, p[3]
        validation_set <- lapply(blocksH0, function(x){
          sample(x$id, size = min(max_size_blockH0_validation, length(x$id)), replace = replace)})
        # sample only the proportion set for validation
        val <- sort(sample(unname(unlist(validation_set)), size = N[3]))
        if(length(val) > N[3]) {
          val <- sort(sample(val, N[3]))
        } else {
          if(length(val) < N[3]) {
            warning(paste0("The size of the validation set is smaller than intended; you should check the values of p[3] (", p[3]
                           ,") and/or max_size_blockH0_validation (", max_size_blockH0_validation, ").
                           As a reference, the number of points assigned to validation are", length(val), ", and the number of points wanted for validation are ", N[3], "."))
          }
        }

        # table(validation_set, spatial$blockH0)

        # set blocks H1 and H2 to be sampled for train/fit and test

        # remove observations already included for validation
        index_train_test <- sp_strat$id |>
          setdiff(val)
        # split by blockH1
        sp_strat_train_test <- sp_strat[sp_strat$id %in% index_train_test,] # remove rows for validation

        # set blocks H1 to be sampled
        blocksH1 <- split(sp_strat_train_test, sp_strat_train_test$blockH1)
        # sample pts in each block H1
        # can we remove the limitation of number of blocks here?
        blocksH1 <- blocksH1[sample.int(length(blocksH1), min(max_number_blocksH1_train, length(blocksH1)))]

        # set blocks H2 within those H1 to be sampled
        # check which blocksH2 are present in the selection H1 blocks above
        all_blocksH2 <- as.numeric(as.character(unlist(lapply(blocksH1, function(x){
          levels(as.factor(x$blockH2)) }))))

        # assign H2 blocks to test/calibration, the rest is for fitting
        # for each H1 block, sample 1 blockH2 for ctest/alibration
        test_blocksH2 <- as.numeric(as.character(unlist(lapply(blocksH1, function(x){
          sample(levels(as.factor(x$blockH2)), 1)})))) # we can change later for more than 1

        # start by setting all H2 blocks in the selected H1 for fitting
        train_set <- sp_strat_train_test[sp_strat_train_test$blockH2 %in% all_blocksH2 &
                                           !(sp_strat_train_test$blockH2 %in% test_blocksH2),]
        # resample data to avoid too many observations of more heavily sampled populations (blockH0 level)
        train_set <- lapply(split(train_set, train_set[[colH0]]), function(x) {
          sample(x$id, size = min(max_size_blockH0_train, length(x$id)), replace = replace) })
        # sample only the desired number of points
        train <- sort(sample(unname(unlist(train_set)), size = N[1]))
        # train/fitting set
        if(length(train) > N[1]) {
          train <- sort(sample(train, N[1]))
        } else {
          if(length(train) < N[1]) {
            train <- c(train, rep(NA, N[1]))[1:N[1]]
            warning("The size of the train/fitting set is smaller than intended; you should check and possibly decrease the value of p[1] and/or increase the value of max_number_blocksH1_train.")
            warning("Replacing the missing test observations by NA. This should be avoided.")
          }
        }

        # test/calibrate
        # set sampled blocksH2 to calibration
        test_set <- sp_strat_train_test[sp_strat_train_test$blockH2 %in% all_blocksH2 & (sp_strat_train_test$blockH2 %in% test_blocksH2),]
        # train/fitting set
        if(nrow(test_set) > N[2]) {
          test <- sort(sample(test_set$id, N[2]))
        } else {
          test <- sort(test_set$id)
          if(length(test) < N[2]) {
            test <- c(test, rep(NA, N[2]))[1:N[2]]
            #print(sum(is.na(test)))
            warning("The size of the test/calibration set is smaller than intended; you should check and possibly decrease the value of p[2] and/or increase the value of max_number_blocksH1_train.")
            warning("Replacing the missing test observations by NA.")
          }
        }

        #subsample fitting data for increased speed
        # fitting_set <- sample(na.omit(spStrat$id[spStrat$set=="fitting"]), size=min(max_fit_size, sum(spStrat$set=="fitting", na.rm=T)))
        # spStrat$set <- ifelse(spStrat$set == "fitting", NA, spStrat$set)
        # spStrat$set <- ifelse(is.na(spStrat$set) & spStrat$id %in% fitting_set, "fitting", spStrat$set)

        # check that no sample is empty
        if(length(train) == 0) stop("No observations were assigned for model fitting.")
        if(length(test) == 0) stop("No observations were assigned for model testing/calibration")
        if(length(val) == 0) stop("No observations were assigned for model validation")

        # update samples
        indexes[[1]][,col] <- train
        indexes[[2]][,col] <- test
        indexes[[3]][,col] <- val
      } # end of loop for times

    } # end of condition for presence of H0
  } # end of spatial stratification condition

  # Organize the output
  # First we had the possibility of returning it in either the list or the matrix type
  # However, we added the blockH0 in the output, so we keep only the list output
  # if (list) {
  #   out <- lapply(indexes, function(x) {
  #     out_samp <- as.data.frame(x, stringsAsFactors = TRUE)
  #     attributes(out_samp) <- NULL
  #     names(out_samp) <- oneimpact::pretty_seq(out_samp)
  #     out_samp
  #   })
  #   names(out) <- c("train", "test", "validate")
  #   out$blockH0 <- colH0
  # }
  # else {
  #   out <- lapply(indexes, function(x) {
  #     out_samp <- x
  #     colnames(out_samp) <- oneimpact::pretty_seq(1:ncol(out_samp))
  #     out_samp
  #   })
  #   names(out) <- c("train", "test", "validate")
  # }
  out <- lapply(indexes, function(x) {
    out_samp <- as.data.frame(x, stringsAsFactors = TRUE)
    attributes(out_samp) <- NULL
    names(out_samp) <- oneimpact::pretty_seq(out_samp)
    out_samp
  })
  names(out) <- c("train", "test", "validate")

  if(length(colH0) == 1) {
    out$blockH0 <- blocks[[colH0]]
  } else {
    out$blockH0 <- colH0
  }

  if(is.null(colH0)) out$blockH0 <- NULL
  if(is.null(sp_strat)) out$sp_strat_id <- NULL else
    out$sp_strat_id <- sp_strat$id

  out
}

#' Create pretty, standardized sequences from numbers
#'
#' @param x `[vector,numeric]` \cr Vector of integers.
#' @param name `[character(1)="Resample"]` \cr String to be
#' added to the numbers in the sequence.
#'
#' @examples
#' pretty_seq(1:20)
#'
#' @keywords internal
#' @export
pretty_seq <- function (x, name = "Resample")
  paste(name, gsub(" ", "0", format(seq(along = x))), sep = "")
