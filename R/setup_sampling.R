#' Setup blocks of points for fitting, calibration, and validation using spatial stratification
#'
#' H0: validation
#' H1, H2
#'
#' @param spStrat A `data.frame` with spatially stratified levels in blocks H0 (population,
#' individual, or study area), H1 (level 1) and H2 (level 2). Each data point is associated
#' to a different set of blocks.
#' @param seed Seed for random sampling.
#' @param sampling_params Sampling parameters, the output from [oneimpact::set_sampling_params].
#'
#' @return A vector with length equals to the number of positions, setting whether the
#' points are set to fitting (training), calibration (testing), and validation.
#'
#' @examples
#' library(terra)
#' library(amt)
#' data(reindeer)
#'
#' # set spatially stratified blocks
#' spst <- spat_strat(reindeer, block_size = 10000, coords = c("x", "y"))
#' # set sampling parameters
#' sampling_params <- set_sampling_params(1000, 10, 1000)
#' # set train, test, and validation sets
#' set_validation_sampling(spst, seed = 1234, sampling_params)
#' # comment
#'
#' @export
set_validation_sampling <- function(spStrat, seed = 1234, sampling_params){

  # Create set column
  spStrat$set <- NA_character_

  # set seed
  set.seed(seed)

  # validation sample from populations (i.e. blockH0)
  # split by pop
  blocksH0 <- split(spStrat, spStrat$blockH0)
  # sample pts in each pop
  validation_set <- unlist(lapply(blocksH0, function(x){
    sample(x$id, min(sampling_params$max_validation_size_blockH0, length(x$id)))}))
  # assign original points to a validation set
  spStrat$set <- ifelse(spStrat$id %in% validation_set, "validation", spStrat$set)

  # sample calibration and fitting for H1 blocks
  # split by blockH1
  blocksH1 <- split(spStrat[is.na(spStrat$set),], spStrat$blockH1[is.na(spStrat$set)])
  # sample pts in each block H1
  blocksH1 <- blocksH1[sample.int(length(blocksH1), min(sampling_params$max_number_blockH1, length(blocksH1)))]
  # check which blocksH2 are present in the selection H1 blocks above
  all_blocksH2 <- as.numeric(as.character(unlist(lapply(blocksH1, function(x){
    levels(as.factor(x$blockH2)) }))))

  # assign H2 blocks to calibration, the rest is for fitting
  # for each H1 block, sample 1 blockH2 for callibration
  calibration_blocksH2 <- as.numeric(as.character(unlist(lapply(blocksH1, function(x){
    sample(levels(as.factor(x$blockH2)), 1)})))) # we can change later for more than 1
  # start by setting all H2 blocks in the selected H1 for fitting
  spStrat$set <- ifelse(is.na(spStrat$set) & spStrat$blockH2 %in% all_blocksH2, "fitting", spStrat$set)
  # set sampled blocksH2 to calibration
  spStrat$set <- ifelse(spStrat$set == "fitting" & spStrat$blockH2 %in% calibration_blocksH2, "calibration", spStrat$set)

  #subsample fitting data for increased speed
  fitting_set <- sample(na.omit(spStrat$id[spStrat$set=="fitting"]), size=min(sampling_params$max_fit_size, sum(spStrat$set=="fitting", na.rm=T)))
  spStrat$set <- ifelse(spStrat$set == "fitting", NA, spStrat$set)
  spStrat$set <- ifelse(is.na(spStrat$set) & spStrat$id %in% fitting_set, "fitting", spStrat$set)

  return(spStrat$set)
}
