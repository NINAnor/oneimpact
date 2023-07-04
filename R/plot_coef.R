#' Plot the coefficients of bags of models
#'
#' Options - one term, multiple bars side by side
#' Colors for negative and positive
#'
#' @export
plot_coef <- function(bag, terms = "all",
                      models = 1:bag$n,
                      weighted = TRUE,
                      what = c("all_models", "average")[1]) {

  # get the coefficients for a subset of models
  coef <- bag$coef[, models]
  w <- bag$weights[models]

  # get the terms of interest
  if(terms != "all") {

    if(is.numeric(terms)) {
      coef <- coef[terms,]
    } else {
      coef <- coef[grepl(terms, rownames(coef)),]
    }

  }

  # weigh (or not) the coefficients
  if(weighted) {

    # all models
    if(what == "all_models") {
      wc <- coef * w
    } else {

      # weighted average
      if(what == "average") {
        wc <- coef %*% w
      }
    }

  } else {
    # raw coefficients, not weighted

    # all models
    if(what == "all_models") {
      wc <- coef
    } else {

      # average
      if(what == "average") {
        ww <- w
        ww[] <- 1/length(w)

        wc <- coef %*% ww
      }
    }

  }

  # data frame and reshape
  if(ncol(coef) > 1) {
    # data frame
    df <- data.frame(x = rownames(coef), wc)

    # reshape
    df <- df |>
      tidyr::pivot_longer(cols = 2:ncol(df), names_to = "resample", values_to = "y")
  } else {
    # data frame
    df <- data.frame(x = rownames(coef), y = wc)
  }

  # plot
  gg <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_bar(stat="identity", fill="steelblue") +
    # scale_y_continuous(trans='log10') +
    ggplot2::theme_minimal() +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Variable", y = "Coefficient")

  if(ncol(coef) > 1) gg <- gg + ggplot2::facet_wrap(~ resample)

  print(gg)
}
