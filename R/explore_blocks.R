#' Explore hierarchical blocks before and after spatial stratification
#'
#' Functions to explore the number of cases and observations for the different
#' sampling units possibly used as the base H0 hierarchical level, such as population ID,
#' study area, animal ID, or year.
#'
#' @param data `[data.frame]` \cr Data.
#' @param colH0 `[character]` \cr Name of the column to be used as the H0 level.
#'
#' @name explore_blocks
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
    dplyr::summarise(n = n(),
                     n_presences = sum(.data[[col_case]] == 1),
                     # n_strata = length(unique(.data[[col_strata]])),
                     .by = all_of(groups))
}

# explore_blocks_pre(dat, colH0 = "studyarea",
#                    animal_id = "animal_id",
#                    col_case = "case_")

# explore_blocks_pre(dat, colH0 = "year",
#                    animal_id = "animal_id",
#                    col_case = "case_")

#' @rdname explore_blocks
#' @export
explore_blocks <- function(blocks) {

  # H0
  blockH0_n <- length(unique(blocks$blockH0))
  blocksH0 <- split(blocks, blocks$blockH0)
  blockH0_size_blocks <- sapply(blocksH0, nrow)

  # H1
  blockH1_n <- length(unique(blocks$blockH1))

  list(blockH0_n = blockH0_n,
       blockH0_size_blocks = blockH0_size_blocks,
       blockH1_n = blockH1_n)
}
# blocksH0 <- split(sp_strat, sp_strat$blockH0)
# sapply(blocksH0, length)
#histogram
# get average and mean number of observations per level
