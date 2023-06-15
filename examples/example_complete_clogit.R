#------------------------
# example of SSF

# packages
library(oneimpact)
library(glmnet)
library(terra)
data("reindeer_annotated")

data_annotated <- tibble::as_tibble(reindeer_annotated)

# plot

# used vs available
data_annotated_v <- terra::vect(data_annotated, geom = c("x2_", "y2_"))
plot(data_annotated_v[data_annotated_v$case_ == F], col = "grey", cex = 0.3)
points(data_annotated_v[data_annotated_v$case_ == T], col = "red", cex = 0.3)

# individuals, used points
plot(data_annotated_v[data_annotated_v$case_ == 1],
     col = data_annotated_v[data_annotated_v$case_ == 1]$animal_year_id, cex = 0.3)

# use add_zoi_formula() and filter_na_strata()
f <- case_ ~ strata(step_id) + sl_*start_roadsXXX + sl_*start_cabinsXXX +
  end_roadsXXX + end_cabinsXXX
f <- add_zoi_formula(f, zoi_radius = c(1000, 2000, 3000), pattern = "XXX")$formula
filter_na_strata(f, data_annotated)

# create blocks
# H0 = individuals, random sampling
data_case1 <- data_annotated |>
  dplyr::filter(case_ == 1)

# whole dataset size
(size_case1 <- length(data_case1$animal_year_id))
# number of blocks H0 (individuals)
(n_blocks <- length(unique(data_case1$animal_year_id)))
# number of observation per block
(obs_block <- table(data_case1$animal_year_id))
# variation in observation per block
stats::quantile(obs_block)

samples <- data_case1 |>
  create_resamples(y = data_case1$case_,
                   p = c(0.2, 0.2, 0.2), times = 10,
                   max_size_validation_blockH0 = 80,
                   max_number_fit_blockH1 = 10,
                   sp_strat = NULL,
                   colH0 = data_case1$animal_year_id)

# check blocks
blocksH0 <- split(samples$validate[[1]], samples$blockH0[samples$validate[[1]]])
sapply(blocksH0, length)

# get variable referring to strata
strat <- extract_response_strata(f)$strata

# fit net_clogit
# not run
if(FALSE) {
  # train and test
  train_data  <- data_annotated[data_annotated[[strat]] %in% data_case1[[strat]][samples$train[[1]]], ]
  test_data <- data_annotated[data_annotated[[strat]] %in% data_case1[[strat]][samples$test[[1]]], ]
  validate_data <- data_annotated[data_annotated[[strat]] %in% data_case1[[strat]][samples$validate[[1]]], ]


  fitted <- net_clogit(f, train_data)
  fitted

  # get variables
  wcols <- extract_response_strata(f, other_vars = TRUE)
  f2 <- as.formula(paste0(wcols$response, " ~ -1 + ", wcols$other_vars))

  #----
  # Variable selection step

  # Calibration with conditionalBoyce index
  # does not work fro cv.glmnet
  # predict.glmnet(fit, test_data, type = "response")??
  pred_vals <- model.matrix(f2, test_data) %*% coef(fitted) # multiple fits?
  # compute conditional Boyce for all lambda parameters
  d <- apply(pred_vals, 2, function(x = x, y = y, strat = strat){
    metric(data.frame(x = x, y = y, strat = strat), errors=F)},
    y = test_data[[wcols$response]], strat = test_data[[wcols$strata]])

  which.max(d)
  plot(fitted$lambda, d, type = "o")
  fitted$lambda[which.max(d)]
  fitted$beta[,which.max(d)]
  data.frame(coef(fitted)[,which.max(d)])

  # initiate results object
  results <- list()
  results$lambda <- fitted$lambda[which.max(d)] # lambda
  results$coef <- matrix(coef(fitted)[,which.max(d)]) # coefficients
  rownames(results$coef) <- names(coef(fitted)[,which.max(d)])
  results$var_names <- names(coef(fitted)[,which.max(d)]) # variable names

  # get predicted values based on the training and testing data
  train_pred_vals <- model.matrix(f2, train_data) %*% results$coef
  test_pred_vals <- model.matrix(f2, test_data) %*% results$coef

  # save results
  results$fit_score <- metric(data.frame(x = train_pred_vals,
                                         y = train_data[[wcols$response]],
                                         strat = train_data[[wcols$strata]]))
  results$calibration_score <- max(d)

  #----
  # Validation step
  val_pred_vals <- model.matrix(f2, validate_data) %*% results$coef
  val <- data.frame(x = val_pred_vals,
                    y = validate_data[[wcols$response]],
                    strat = validate_data[[wcols$strata]])
  # val <- split(val, samples$blockH0[sort(match(val$strat, spStrat$id))])
  val <- split(val, samples$blockH0[match(val$strat, validate_data[[wcols$strata]])])
  results$validation_score <- unlist(lapply(val, metric))

  # Validation habitat only
  pred_vals_kernel <- kernel_prediction(f, validate_data,
                                        kernel_vars = kernel_vars,
                                        coefs = results$coef[,1])

  pred_vals_habitat <- val_pred_vals - pred_vals_kernel # does it make sense??
  hab <- data.frame(x = pred_vals_habitat,
                    y = validate_data[[wcols$response]],
                    strat = validate_data[[wcols$strata]])
  hab <- split(hab, samples$blockH0[match(hab$strat, validate_data[[wcols$strata]])])
  results$habitat_validation_score <- unlist(lapply(hab, metric))
  #plot(results$validation_score, results$habitat_validation_score)
}

# use fit_net_clogit
# not run
if(FALSE) {
  fitted1 <- fit_net_clogit(f, data_annotated,
                            samples = samples, i = 5,
                            kernel_vars = "sl_",
                            metric = AUC)
  fitted1
}
