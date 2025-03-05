#' Remove missing values and ensure case-control per stratum
#'
#' Remove strata for which there are only 1s (used points, presence) or only 0s
#' (available points, background or absence).
#'
#' @param f `[formula]` \cr Formula for model fitting.
#' @param data `[data.frame]` \cr Data set to be checked.
#'
#' @return Cleaned `data.frame`, removing from the input `data` the rows with NA
#' in any of the columns and all the strata for which there are only presences or
#' only absences.
#'
#' @example examples/filter_na_strata_example.R
#'
#' @references [amt::remove_incomplete_strata()]
#'
#' @export
filter_na_strata <- function(f, data){

  # removes missing data in any columns
  # this removes too much data!! for now we do not do it here.
  data_clean <- data #na.omit(data[,all.vars(f)])
  # get strata variable
  strat <- extract_response_strata(f)$strata
  # get case - response variable
  cas <- extract_response_strata(f)$response
  # split data.frame by stratum
  lp <- tidyr::nest(data_clean[,c(cas, strat)], .by = tidyselect::all_of(strat))

  # tic(); lp <- split(data_clean[,c(cas, strat)], data_clean[,strat]); toc()
  # tic(); lp <- tidyr::nest(data_clean[,c(cas, strat)], .by = tidyselect::all_of(strat)); toc() # much faster

  # lp <- split(data_clean[,c(cas, strat)], data_clean[,strat])
  # test <- unlist(lapply(lp, function(x){mean(x[[cas]])}))
  # keep <- unlist(lapply(lp, function(x){x[1,2]}))

  test <- unlist(lapply(lp$data, function(x){mean(x[[cas]])}))
  keep <- lp[[strat]]
  keep <- keep[(test>0 & test<1)]
  data_clean <- data_clean[data_clean[[strat]] %in% keep,]
  return(data_clean)
}
