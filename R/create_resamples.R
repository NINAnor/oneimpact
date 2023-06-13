#' Create samples for fitting, calibrating, and validating models
#'
#' The function creates samples of the data to be included for three different
#' model blocks: train (fitting), validating, and testing models.
#' By default, samples are created through bootstrapping, i.e. with replacement. This means
#' data observations can be repeated within a given sample block, but observations included
#' in one block are necessarily excluded from the other blocks (e.g. observations selected
#' for test will be absent from train and validation blocks).
#'
#' @param y a vector of outcomes.
#' @param times the number of partitions or samples to create
#' @param p a 3 element vector with the percentage of data that goes to training, validation, and test
#' @param list logical - should the results be in lists (TRUE) or matrices with the number of rows equal to
#' floor(p * length(y)) and times columns.
#' @param replace whether to perform the bootstrap sampling with or without replacement (Default = TRUE).
#'
#' @examples
#' y <- runif(200)
#' samples <- create_resamples(y, p = c(0.4, 0.2, 0.2), times = 5)
#' samples
#'
#' @export
create_resamples <- function (y, times = 10,
                              p = c(0.4, 0.2, 0.2),
                              list = TRUE,
                              replace = TRUE)
{
  if (inherits(y, "Surv"))
    y <- y[, "time"]

  # set number of points in each
  # here we can use different rules
  N <- floor(p * length(y))

  # set matrices
  # 3 matrices, for train, test, and validate
  indexes <- lapply(N, function(n) matrix(0, ncol = times, nrow = n))

  # sample
  # first validate, then train and test, mutually exclusive
  for(col in seq(times)) {

    # validation
    index_val <- seq(along = y)
    val <- sort(sample(index_val, size = nrow(indexes[[3]]), replace = TRUE))

    # fitting
    index_fit <- index_val |>
      setdiff(val)
    fit <- sort(sample(index_fit, size = nrow(indexes[[1]]), replace = TRUE))

    # calibrating
    index_cal <- index_fit |>
      setdiff(fit)
    cal <- sort(sample(index_cal, size = nrow(indexes[[2]]), replace = TRUE))

    indexes[[1]][,col] <- fit
    indexes[[2]][,col] <- cal
    indexes[[3]][,col] <- val
  }

  if (list) {
    out <- lapply(indexes, function(x) {
      out_samp <- as.data.frame(x, stringsAsFactors = TRUE)
      attributes(out_samp) <- NULL
      names(out_samp) <- pretty_seq(out_samp)
      out_samp
    })
    names(out) <- c("train", "validate", "test")
  }
  else {
    out <- lapply(indexes, function(x) {
      out_samp <- x
      colnames(out_samp) <- pretty_seq(1:ncol(out_samp))
      out_samp
    })
    names(out) <- c("train", "validate", "test")
  }
  out
}

# Helper function to create resample names
pretty_seq <- function (x, name = "Resample")
  paste(name, gsub(" ", "0", format(seq(along = x))), sep = "")
