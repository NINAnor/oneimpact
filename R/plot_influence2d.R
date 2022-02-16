#' Plot nearest and cumulative influence in 2 dimensions
#'
#' To add: options for normalization, prob etc.
#'
#' @param points Vector of values in the x axis where the infrastructure (or the origin
#' of the functions) are located
#'
#' @return a ggplot object with the nearest or cumulative influence plot.
#' If `return_df = TRUE`, it returns a list with the ggplot object and the `data.frame`
#' with values.
#'
#' @example examples/plot_influence2d_example.R
#'
#' @export
plot_influence2d <- function(points,
                             zoi,
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
  dat_y <- purrr::map_dfc(points, ~ fun(x = x, zoi = zoi, origin = .x, oneside = FALSE))

  # should functions accumulate or not?
  if(cumulative) y <- apply(dat_y, 1, sum, na.rm = na.rm) else
    y <- apply(dat_y, 1, max, na.rm = na.rm)

  # tibble
  dat <- tibble::tibble(x = x, y = y)

  # plot
  g1 <- ggplot(dat, aes(x = x, y = y)) +
    geom_line() +
    labs(x = "Distance", y = "Influence") +
    theme_bw()

  # return plot (and maybe df)
  if (isTRUE(return_df)) {
    return(list(influence_plot = g1, influence_df = dat))
  } else {
    return(g1)
  }

}
