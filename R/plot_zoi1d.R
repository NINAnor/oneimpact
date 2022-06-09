#' Plot nearest and cumulative zone of influence function in 1 dimension
#'
#' Plots the zone of influence (ZoI) functions in 1 dimensional space, for illustration.
#' When there is more than one value for `points` (the location of infrastructure
#' or sources of disturbance), either the ZoI of the nearest feature or the
#' cumulative ZoI can be plotted. The ZoI of the nearest feature corresponds to the maximum ZoI value from all
#' infrastructure at each position. The cumulative ZoI corresponds to the sum of the ZoI
#' of all infrastructure at each position.
#'
#' In practice, the `plot_zoi1d` computes the ZoI value for each feature
#' whose locations in 1 dimension are defined by `points` and calculates the maximum
#' (ZoI of the nearest) or sum (cumulative ZoI) of all values This is done
#' for a series of points in 1 dimensional space in the range `range_plot` (with
#' steps defined by `step`) and finally plotted.
#'
#' To add: options for normalization, prob etc.
#'
#' @param points `[numeric]` \cr Vector of values in the x axis where the infrastructure are located,
#' in 1 dimension
#' (more broadly, the location of sources of disturbance or spatial variables, or the point of
#' origin of the ZoI functions).
#' @param zoi_radius `[numeric(1)]` \cr Zone of influence (ZoI) radius,
#' the distance at which the ZoI vanishes or goes below a given minimum threshold value.
#' See [zoi_functions] for details.
#' @param fun `[function]` \cr A decay function that represents the Zone of Influence (ZoI).
#' Different functions might represent different shapes for the decay of the ZoI.
#' See [zoi_functions] for some examples.
#' @param cumulative `[logical(1)=FALSE]` \cr If `TRUE`, the cumulative ZoI is plotted. If `FALSE`
#' (default), the ZoI of the nearest feature is plotted.
#' @param range_plot `[numeric(2)=c(0,12)]` \cr A vector c(xmin,xmax) with the x range of
#' the ZoI plot.
#' @param step `[numeric(1)=0.01]` \cr Size of the step increment used to define the series
#' of x positions for which the ZoI is computed, within the x range defined by `range_plot`.
#' @param return_df `[logical(1)=FALSE]` \cr If TRUE, a `data.frame` with `x` values and their
#' corresponding ZoI values is returned, besides the plot.
#' @param ... Additional parameters passed to the ZoI decay function `fun`.
#'
#' @return A `ggplot` object with the nearest or cumulative influence plot.
#' If `return_df = TRUE`, it returns a list with the `ggplot` object and the `data.frame`
#' with values for `x` (position in 1d space) and `y` (ZoI value).
#'
#' @example examples/plot_zoi1d_example.R
#'
#' @export
plot_zoi1d <- function(points,
                       zoi_radius,
                       fun = exp_decay,
                       cumulative = FALSE,
                       range_plot = c(0, 12),
                       step = 0.01,
                       na.rm = TRUE,
                       return_df = FALSE,
                       ...) {

  # if function is a character
  if(is.character(fun)) {
    fun <- get(fun)
  }

  # create x
  x <- seq(range_plot[1], range_plot[2], step)
  # apply function to each point
  dat_y <- purrr::map_dfc(points, function(z) fun(x = x, zoi_radius = zoi_radius, origin = z, oneside = FALSE, ...))

  # should functions accumulate or not?
  if(cumulative) y <- lapply(dat_y, sum, na.rm = na.rm) else
    y <- apply(dat_y, 1, max, na.rm = na.rm)

  # tibble
  dat <- tibble::tibble(x = x, y = y)

  # plot
  g1 <- ggplot(dat, aes(x = x, y = y)) +
    geom_line() +
    labs(x = "Distance", y = "Zone of Influence") +
    theme_bw()

  # return plot (and maybe df)
  if (isTRUE(return_df)) {
    return(list(influence_plot = g1, influence_df = dat))
  } else {
    return(g1)
  }

}
