#' Plot the functions for the nearest and cumulative zone of influence in 1 dimension
#'
#' Plots the functions to represent the zone of influence (ZOI) in 1 dimensional space, for illustration.
#' When there is more than one value for `points` (the location of infrastructure
#' or sources of disturbance), two metrics might be plotted: the ZOI of the nearest feature or the
#' cumulative ZOI. The ZOI of the nearest feature corresponds to the maximum ZOI value from all
#' infrastructure at each position. The cumulative ZOI corresponds to the sum of the ZOI
#' of all infrastructure at each position.
#'
#' In practice, the `plot_zoi1d()` computes the ZOI value for each feature
#' whose locations in 1 dimension are defined by `points` and calculates the maximum
#' (ZOI of the nearest) or sum (cumulative ZOI) of all values. This is done
#' for a series of points in 1 dimensional space in the range `range_plot` (with
#' steps defined by `step`) and plotted.
#'
#' To add: options for normalization, prob etc.
#'
#' @param points `[numeric]` \cr Vector of values in the x axis where the infrastructure are located,
#' in 1 dimension
#' (more broadly, the location of sources of disturbance or spatial variables, or the point of
#' origin of the ZOI functions).
#' @param radius `[numeric(1)]` \cr Radius of the zone of influence (ZOI),
#' the distance at which the ZOI vanishes or goes below a given minimum limit value.
#' See [zoi_functions] for details.
#' @param fun `[function]` \cr A decay function that represents the Zone of Influence (ZOI).
#' Different functions might represent different shapes for the decay of the ZOI.
#' See [zoi_functions] for some examples.
#' @param zoi_metric `[character(1)="nearest"]{"nearest", "cumulative"}` \cr
#' Which metric of zone of influence should be plotted. If `"nearest"`(default), the ZOI of
#' the nearest feature is plotted. If `"cumulative"`, the cumulative ZOI is plotted.
#' @param range_plot `[numeric(2)=c(0,12)]` \cr A vector `c(xmin, xmax)` with the x range of
#' the ZOI plot.
#' @param step `[numeric(1)=0.01]` \cr Size of the step increment used to define the series
#' of x positions for which the ZOI is computed, within the x range defined by `range_plot`.
#' @param return_df `[logical(1)=FALSE]` \cr If TRUE, a `data.frame` with `x` values and their
#' corresponding ZOI values is returned, besides the plot.
#' @param na.rm `[logical(1)=TRUE]` \cr Should `NA` values be removed when computing
#' either the sum (for the cumulative ZOI) or the maximum (for the ZOI of the nearest
#' feature)?
#' @param ... Additional parameters passed to the ZOI decay function `fun`.
#'
#' @return A `ggplot` object with the nearest or cumulative influence plot.
#' If `return_df = TRUE`, it returns a list with the `ggplot` object and a `data.frame`  with
#' values for `x` (position in 1d space) and `y` (ZOI value).
#'
#' @example examples/plot_zoi1d_example.R
#'
#' @export
plot_zoi1d <- function(points,
                       radius,
                       fun = exp_decay,
                       zoi_metric = c("nearest", "cumulative")[1],
                       range_plot = c(0, 12),
                       step = 0.01,
                       na.rm = TRUE,
                       return_df = FALSE,
                       ...) {

  # if function is a character
  if(is.character(fun)) {
    fun <- get(fun, asNamespace("oneimpact"))
  }

  # create x
  x <- seq(range_plot[1], range_plot[2], step)
  # apply function to each point
  dat_y <- sapply(points, function(z) fun(x = x, radius = radius, origin = z, oneside = FALSE, ...)) |>
    as.data.frame()

  # should functions accumulate or not?
  if(zoi_metric == "cumulative") y <- apply(dat_y, 1, sum, na.rm = na.rm) else
    y <- apply(dat_y, 1, max, na.rm = na.rm)

  # tibble
  dat <- tibble::tibble(x = x, y = y)

  # plot
  g1 <- ggplot2::ggplot(dat, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Distance", y = "Zone of Influence") +
    ggplot2::theme_bw()

  # return plot (and maybe df)
  if (isTRUE(return_df)) {
    return(list(zoi_plot = g1, zoi_df = dat))
  } else {
    return(suppressWarnings(g1))
  }

}
