#' Explore potential hierarchical blocks before sampling or spatial stratification
#'
#' Function to explore the number of cases and observations for the different
#' sampling units possibly used as the base H0 hierarchical level, such as population ID,
#' study area, animal ID, or year, before spatial stratification or creating samples
#' for the bootstrapped approach.
#' The function can help understand how imbalanced is the data across H0 levels used
#' for validation.
#'
#' @param data `[data.frame,tibble]` \cr Complete data set to be analyzed.
#' @param colH0 `[character]` \cr Name of the column in `data` to be used as the H0
#' hierarchical level, intended for model validation.
#' @param animal_id `[character]` \cr Name of the column in `data` representing
#' animal ID. If `NULL` (default), summaries are not created for individuals.
#' @param col_case `[string(1)="case"]` \cr Name of the column in `data` representing
#' the case or used/available points. Default is `"case"`.
#'
#' @example examples/explore_blocks_pre_example.R
#'
#' @name explore_blocks_pre
#' @export
explore_blocks_pre <- function(data,
                               colH0,
                               animal_id = NULL,
                               col_case = "case") {

  # groups
  if(is.null(animal_id)) {
    groups <- c(colH0)
  } else {
    groups <- c(colH0, animal_id)
  }

  # summary
  data |>
    dplyr::summarise(n = dplyr::n(),
                     n_presences = sum(.data[[col_case]] == 1),
                     # n_strata = length(unique(.data[[col_strata]])),
                     .by = all_of(groups))
}

#' Explore hierarchical blocks after sampling or spatial stratification
#'
#' Function to explore the number of cases and observations for the
#' sampling units used as the base H0 hierarchical level used for model validation,
#' such as population ID, study area, animal ID, or year, after spatial stratification.
#' It also allows exploring the the spatially stratified levels H1 and H0, to be used
#' for model fitting and tuning (train and test), respectively.
#' The function can help understand how imbalanced are data across H0, H1, and H2
#' levels used for validation, training/fitting, and tuning/testing, as well as how many
#' blocks exist there are for each type and their sizes.
#'
#' @param blocks `[data.frame]` \cr A `data.frame` returned by the function
#' [oneimpact::spat_strat()], which a list of blocks H0 used for validation and
#' blocks H1 and H2 spatially stratified.
#'
#' @returns A list with:
#' - blockH0_n: the number of levels/blocks in the hierarchical level H0, to be
#' used for validation;
#' - blockH0_size_blocks: the number of (used) observations in each block of
#' the hierarchical level H0, to be used for validation;
#' - blockH1_n: the number of different blocks (spatial strata) in the hierarchical
#' level H1, to be used for model fitting.
#' - blockH1_size_blocks: the number of (used) observations in each block of
#' the hierarchical level H1, to be used for model fitting;
#' - blockH2_n: the number of different blocks (spatial strata) in the hierarchical
#' level H2, to be used for model tuning.
#' - blockH2_size_blocks: the number of (used) observations in each block of
#' the hierarchical level H2, to be used for model tuning.
#'
#' @example examples/explore_blocks_example.R
#'
#' @name explore_blocks
#' @export
explore_blocks <- function(blocks) {

  # H0
  blockH0_n <- length(unique(blocks$blockH0))
  blocksH0 <- split(blocks, blocks$blockH0)
  # Size of blocks H0
  blockH0_size_blocks <- sapply(blocksH0, nrow)

  # H1
  blockH1_n <- length(unique(blocks$blockH1))
  blocksH1 <- split(blocks, blocks$blockH1)
  # Size of blocks H1
  blockH1_size_blocks <- sapply(blocksH1, nrow)

  # H2
  blockH2_n <- length(unique(blocks$blockH2))
  blocksH2 <- split(blocks, blocks$blockH2)
  # Size of blocks H2
  blockH2_size_blocks <- sapply(blocksH2, nrow)

  list(blockH0_n = blockH0_n,
       blockH0_size_blocks = blockH0_size_blocks,
       blockH1_n = blockH1_n,
       blockH1_size_blocks = blockH1_size_blocks,
       blockH2_n = blockH2_n,
       blockH2_size_blocks = blockH2_size_blocks)
}
# blocksH0 <- split(sp_strat, sp_strat$blockH0)
# sapply(blocksH0, length)
#histogram
# get average and mean number of observations per level
