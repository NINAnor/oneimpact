#' Plot the coefficients of bags of models
#'
#' Options - one term, multiple bars side by side
#' Colors for negative and positive
#'
#' @export
plot_coef <- function(bag, terms = "all",
                      models = 1:bag$n,
                      weighted = TRUE,
                      what = c("all_models", "average")[1],
                      plot_type = c("bars", "points", "histogram")[1],
                      remove_low = 0, remove_high = Inf,
                      order_zoi_radius = FALSE,
                      show_legend = FALSE) {

  # get the coefficients for a subset of models
  coef <- bag$coef[, models]
  w <- bag$weights[models]

  # get the terms of interest
  if(terms != "all") {

    if(is.numeric(terms)) {
      coef <- coef[terms,]
    } else {
      coef <- coef[grepl(terms, rownames(coef)),]

      # check if it is ZOI variables and they should be ordered
      if(order_zoi_radius) {
        radii <- as.numeric(gsub("\\D", "", rownames(coef)))
        coef <- coef[order(radii),]
      }
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

        wc_quant <- sapply(1:nrow(coef), function(i) DescTools::Quantile(coef[i,], weights = w, probs = c(0.025, 0.975), type = 5))
        wc_min <- wc_quant[1,]
        wc_max <- wc_quant[2,]
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
        ww[] <- 1/length(ww)

        wc <- coef %*% ww
        wc_quant <- sapply(1:nrow(coef), function(i) quantile(coef[i,], probs = c(0.025, 0.975), na.rm = TRUE))
        wc_min <- wc_quant[1,]
        wc_max <- wc_quant[2,]
      }
    }

  }

  # data frame and reshape
  if(ncol(coef) > 1 & what != "average") {
    # data frame
    df <- data.frame(x = rownames(coef), wc)

    # reshape
    df <- df |>
      tidyr::pivot_longer(cols = 2:ncol(df), names_to = "resample", values_to = "y")
  } else {
    if(what == "average") {
      # data frame
      df <- data.frame(x = rownames(coef), y = as.numeric(wc[,1]), lower = wc_min, upper = wc_max)
    } else {
      # data frame
      df <- data.frame(x = rownames(coef), y = as.numeric(wc[,1]))
    }
  }

  # negative, positive, zero
  df$signal <- ifelse(df$y > 0, "positive", ifelse(df$y < 0, "negative", "null"))
  # filter thresholds
  df <- df[abs(df$y) >= remove_low & abs(df$y) < remove_high,]
  df$x <- factor(df$x, levels = unique(df$x))

  # plot type
  if(plot_type == "bars") {
    plot_func <- function(df, ...) {
      ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = signal)) +
        ggplot2::geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
        ggplot2::geom_bar(stat = "identity", ...) +
        ggplot2::labs(x = "Variable", y = "Coefficient", fill = "Signal") +
        ggplot2::coord_flip()
    }
  } else {

    if(plot_type == "points") {
      if(what == "all_models") {
        plot_func <- function(df, ...) {
          ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
            ggplot2::geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
            ggplot2::geom_segment(ggplot2::aes(x = x, xend = x, y = 0, yend = y), color = "grey", ...) +
            ggplot2::geom_point(aes(color = signal), size=2) +
            ggplot2::labs(x = "Variable", y = "Coefficient", color = "Signal") +
            ggplot2::coord_flip()
        }
      } else {
        plot_func <- function(df, ...) {
          ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
            ggplot2::geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
            ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.2, color = "grey", ...) +
            ggplot2::geom_point(aes(color = signal), size=2) +
            ggplot2::labs(x = "Variable", y = "Coefficient", color = "Signal") +
            ggplot2::coord_flip()
        }
      }
    } else {

      if(plot_type == "histogram") {
        plot_func <- function(df, ...) {
          ggplot2::ggplot(df, ggplot2::aes(x = y, fill = signal)) +
            ggplot2::geom_vline(xintercept = 0, color = "grey", linetype = "dashed") +
            ggplot2::geom_histogram(position = position_dodge()) +
            ggplot2::labs(x = "Coefficient", y = "Number of models", fill = "Signal")
        }
        if(what == "average") warning("A histogram based on a single average model is not very meaningful.")
      }

    }

  }

  # plot
  gg <- plot_func(df) +
    ggplot2::theme_minimal()

  # options
  # legend
  if(!show_legend) gg <- gg + theme(legend.position = "none")

  # facets per model
  if(ncol(coef) > 1 & plot_type != "histogram" & what == "all_models") {
    gg <- gg + ggplot2::facet_wrap(~resample)
  }

  # facets per variable
  if(nrow(coef) > 1 & plot_type == "histogram") {
    gg <- gg + ggplot2::facet_wrap(~x)
  }

  print(gg)
}

#'
plot_weights <- function(x, pattern = "*", remove_low = 0, remove_high = Inf, normalize = FALSE) {

  # weighted coefs
  w_coef <- x$coef %*% x$weights

  # subset
  w_coef <- w_coef[grepl(pattern, rownames(w_coef)),]

  # normalization
  if(normalize)
    wgt_coef <- w_coef/max(w_coef) else
      wgt_coef <- w_coef

  # data frame
  df <- data.frame(var = factor(names(w_coef), levels = names(wgt_coef), ordered = TRUE),
                   coef = wgt_coef)
  # filter thresholds
  df <- df[abs(df$coef) >= remove_low & abs(df$coef) < remove_high,]

  # plot
  p <- ggplot2::ggplot(data = df, ggplot2::aes(x = var, y = coef)) +
    ggplot2::geom_bar(stat="identity", fill="steelblue") +
    # scale_y_continuous(trans='log10') +
    ggplot2::theme_minimal() +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Variable", y = "Weighted coefficients")
  print(p)
}
