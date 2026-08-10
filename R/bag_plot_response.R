#' Plot response curves from a bag of models
#'
#' This function predicts and plots response curves from a bag of models while
#' varying one or more focal predictors in `dfvar`. Non-focal predictors are kept
#' at baseline values. It supports confidence intervals from weighted quantiles
#' and optional individual-model predictions.
#'
#' The function `plot_response` uses the `predict` to produce the predictions.
#'
#' @param x `[bag,list]` \cr A bag of models created with [oneimpact::bag_models()].
#' @param dfvar `[data.frame]` \cr Data frame with values of focal predictors to vary
#' and predict over. The column names of the `data.frame` might correspond exactly to the model covariates or to
#' parts of that (for instance, "roads_paved_" to refer to all ZOI variables related to paved roads).
#' @param data `[data.frame]` \cr Original data used for model fitting. Used only for
#' taking the categories of the categorical variables. Irrelevant if there is no categorical variables.
#' @param type `[character(1)="linear"]{"linear", "exponential", "logit", "cloglog"}` \cr
#' Prediction scale.
#' @param zoi_shape `[character(1)="exp_decay"]{"exp_decay", "gaussian_decay", "linear_decay", "threshold_decay"}` \cr
#' ZOI decay shape used when `zoi = TRUE`.
#' @param which_cumulative `[character(1)="cumulative"]` \cr Pattern used to identify
#' cumulative ZOI terms.
#' @param ci `[logical(1)=TRUE]` \cr If `TRUE` (default), plot weighted confidence intervals
#' using `wq_probs`.
#' @param indiv_pred `[logical(1)=FALSE]` \cr If `TRUE`, include curves from individual
#' models (only models with positive weight).
#' @param wq_probs `[numeric,vector=c(0.025, 0.5, 0.975)]` \cr Weighted quantile
#' probabilities used for confidence intervals and median.
#' @param baseline `[character(1)="median"]{"median", "mean", "zero"}` \cr Baseline
#' value strategy for non-focal predictors. Variable are either kept constant at the mean
#' or median values, or left as zero. Categorical variables are set to their reference
#' level, retrieved from `data`.
#' @param zoi `[logical(1)=FALSE]` \cr If `TRUE`, variables values in `dfvar` are
#' interpreted as distance from a disturbance source, to be transformed into ZOI predictors.
#' @param type_feature_recompute `[logical(1)=FALSE]` \cr If `TRUE`, recompute
#' line- or area-feature raster representation for ZOI calculations.
#' @param type_feature `[character(1)="point"]{"point", "line", "area"}` \cr Feature type
#' for ZOI prediction.
#' @param zoi_limit `[numeric(1)=0.05]` \cr Lower influence threshold used by non-vanishing
#' ZOI functions. See [oneimpact::zoi_functions()].
#' @param resolution `[numeric(1)=100]` \cr Raster resolution used for line-feature ZOI
#' approximation. Used when recomputing the ZOI variables for line and area features.
#' @param line_value `[numeric(1)=1]` \cr Value assigned to line raster cells when
#' `type_feature = "line"`. Used when recomputing the ZOI variables for line and area features.
#' @param ggplot `[logical(1)=TRUE]` \cr If `TRUE`, return a ggplot object; otherwise
#' return a prediction data frame.
#' @param plot_mean `[logical(1)=TRUE]` \cr Plot weighted mean response line.
#' @param plot_median `[logical(1)=TRUE]` \cr Plot weighted median response line.
#' @param n_features `[numeric(1)=1]` \cr Number of features used in ZOI prediction.
#' To represent the cumulative impact of multiple features, use `n_features > 1`.
#' @param normalize `[logical or character]` \cr Optional y-axis normalization:
#' one of `FALSE`, `"mean"`, `"median"`, or `"ci"`. Default is `FALSE`, which assumes
#' no normalization.
#' @param logx `[logical(1)=FALSE]` \cr If `TRUE`, use log10 scaling on x-axis.
#' @param ylim `[NULL or ggplot2 scale/coord limits]` \cr Optional y-axis limits.
#' @param y_lab `[character(1)="Relative Selection Strength"]` \cr Y-axis label.
#' @param col_ci `[character(1)="grey"]` \cr Fill color for confidence ribbon.
#' @param col_indiv `[character(1)="grey"]` \cr Color for individual model lines.
#' @param col_mean `[character(1)="black"]` \cr Color for weighted mean line.
#' @param col_median `[character(1)="red"]` \cr Color for weighted median line.
#' @param linewidth_indiv `[numeric(1)=1.2]` \cr Line width for individual model lines.
#' @param linewidth_mean `[numeric(1)=1.2]` \cr Line width for weighted mean line.
#' @param linewidth_median `[numeric(1)=1.2]` \cr Line width for weighted median line.
#' @param alpha_ci `[numeric(1)=0.5]` \cr Alpha transparency for confidence ribbon.
#' @param alpha_indiv `[numeric(1)=0.3]` \cr Alpha transparency for individual model lines.
#'
#' @return If `ggplot = TRUE`, a ggplot object with response curves.
#' If `ggplot = FALSE`, a data frame with `dfvar`, summary predictions, and
#' optionally individual-model predictions when `indiv_pred = TRUE`.
#'
#' @seealso [oneimpact::predict()]
#'
#' @example examples/bag_plot_response_example.R
#'
#' @name plot_response
#' @export
plot_response <- function(x,
                          dfvar,
                          data,
                          type = c("linear", "exponential", "logit", "cloglog")[1],
                          zoi_shape = c("exp_decay", "gaussian_decay", "linear_decay", "threshold_decay")[1],
                          which_cumulative = "cumulative",
                          ci = TRUE,
                          indiv_pred = FALSE,
                          wq_probs = c(0.025, 0.5, 0.975),
                          baseline = c("median", "mean", "zero")[1],
                          zoi = FALSE,
                          #zoi_vals = c(100, 250, 500, 1000, 2500, 5000, 10000),
                          type_feature = c("point", "line", "area")[1],
                          type_feature_recompute = FALSE,
                          zoi_limit = 0.05,
                          resolution = 100, # resolution for the raster created in create_line_feature_zoi
                          line_value = 1, # value set to the linear infrastructure raster created in create_line_feature_zoi
                          ggplot = T,
                          plot_mean = TRUE,
                          plot_median = TRUE,
                          n_features = 1,
                          normalize = c(FALSE, "mean", "median", "ci")[1],
                          logx = FALSE,
                          ylim = NULL,
                          y_lab = "Relative Selection Strength",
                          col_ci = "grey",
                          col_indiv = "grey",
                          col_mean = "black",
                          col_median = "red",
                          linewidth_indiv = 1.2,
                          linewidth_mean = 1.2,
                          linewidth_median = 1.2,
                          alpha_ci = 0.5,
                          alpha_indiv = 0.3) {

  UseMethod("plot_response")

}


#' @rdname plot_response
#' @export
plot_response.bag <- function(x,
                              dfvar,
                              data,
                              type = c("linear", "exponential", "logit", "cloglog")[1],
                              zoi_shape = c("exp_decay", "gaussian_decay", "linear_decay", "threshold_decay")[1],
                              which_cumulative = "cumulative",
                              ci = TRUE,
                              indiv_pred = FALSE,
                              wq_probs = c(0.025, 0.5, 0.975),
                              baseline = c("median", "mean", "zero")[1],
                              zoi = FALSE,
                              # zoi_vals = c(100, 250, 500, 1000, 2500, 5000, 10000),
                              type_feature = c("point", "line", "area")[1],
                              type_feature_recompute = FALSE,
                              zoi_limit = 0.05,
                              resolution = 100,
                              line_value = 1,
                              ggplot = T,
                              plot_mean = TRUE,
                              plot_median = TRUE,
                              n_features = 1,
                              normalize = c(FALSE, "mean", "median", "ci")[1],
                              logx = FALSE,
                              ylim = NULL,
                              y_lab = "Relative Selection Strength",
                              col_ci = "grey",
                              col_indiv = "grey",
                              col_mean = "black",
                              col_median = "red",
                              linewidth_indiv = 1.2,
                              linewidth_mean = 1.2,
                              linewidth_median = 1.2,
                              alpha_ci = 0.5,
                              alpha_indiv = 0.3){

  # predict for new data set

  # predict only for those variables/coefs of interest for ZOI variables
  if(zoi) include <- colnames(dfvar) else include = "all"

  # predict
  pred <- oneimpact::predict(x,
                             newdata = dfvar,
                             data = data,
                             type = type,
                             wmean = TRUE,
                             wq_probs = wq_probs,
                             include = include,
                             baseline = baseline,
                             zoi = zoi,
                             zoi_shape = zoi_shape,
                             which_cumulative = which_cumulative,
                             type_feature = type_feature,
                             type_feature_recompute = type_feature_recompute,
                             n_features = n_features,
                             zoi_limit = zoi_limit,
                             resolution = resolution,
                             line_value = line_value)
  names(pred) <- if(is.null(wq_probs)) c("mean") else c("lower", "mid", "higher", "mean")


  # predict for individual models
  if(indiv_pred) {
    pred_indiv <- oneimpact::predict(x,
                                     newdata = dfvar,
                                     data = data,
                                     type = type,
                                     wmean = FALSE,
                                     wq_probs = NULL,
                                     include = include,
                                     baseline = baseline,
                                     zoi = zoi,
                                     zoi_shape = zoi_shape,
                                     which_cumulative = which_cumulative,
                                     type_feature = type_feature,
                                     type_feature_recompute = type_feature_recompute,
                                     n_features = n_features,
                                     zoi_limit = zoi_limit,
                                     resolution = resolution,
                                     line_value = line_value)[, x$weights > 0]
  }


  if (ggplot){
    if (zoi){
      if (ncol(dfvar) == 1){

        # data for plotting
        df <- data.frame(x = dfvar[,1], y = pred$mid,
                         y_lower = pred$lower, y_upper = pred$higher, y_mean = pred$mean)
        # data for plotting individual models
        if(indiv_pred) df_indiv <- tidyr::gather(data.frame(x = dfvar[,1], pred_indiv),
                                                 key = "model", value = "y", -x)

        # normalize y axis?
        if(normalize != FALSE) {
          if(normalize == "mean") range_y <- range(df$y_mean[-1]) else
            if(normalize == "median") range_y <- range(df$y[-1]) else
              range_y <- range(df[-1,2:ncol(df)])
            df[2:ncol(df)] <- do.call("cbind", lapply(2:ncol(df), function(i) (df[,i] - range_y[1])/(diff(range_y))))
            if(indiv_pred) df_indiv$y <- (df_indiv$y - range_y[1])/(diff(range_y))
        }

        # plot
        plt <- ggplot2::ggplot(df)

        # confidence interval
        if(ci) {
          plt <- plt + ggplot2::geom_ribbon(ggplot2::aes(x = x,
                                                         ymin = y_lower,
                                                         ymax = y_upper),
                                            fill = col_ci,
                                            alpha = alpha_ci)
        } else {
          # individual lines
          if(indiv_pred) {
            plt <- plt + ggplot2::geom_line(ggplot2::aes(x = x,
                                                         y = y,
                                                         group = model),
                                            data = df_indiv,
                                            color = col_indiv,
                                            linewidth = linewidth_indiv,
                                            alpha = alpha_indiv)
          }
        }
        # median
        if(plot_median)
          plt <- plt + ggplot2::geom_line(ggplot2::aes(x=x,
                                                       y = y),
                                          color = col_mean,
                                          linewidth = linewidth_mean)
        # mean
        if(plot_mean)
          plt <- plt + ggplot2::geom_line(ggplot2::aes(x=x,
                                                       y = y_mean),
                                          color = col_median,
                                          linewidth = linewidth_median)

        plt <- plt + ggplot2::labs(x = ifelse(zoi, "Distance (m)", names(dfvar)),
                                   y = y_lab, title = "")
      }
      # if (ncol(dfvar)==2){
      #   df <- data.frame(x=dfvar[,1], grp=as.factor(newdata[,names(dfvar)[2]]), y=pred$mid, y_lower = pred$lower, y_upper=pred$higher, y2=pred$mean)
      #   plt <- ggplot(df, aes(x=x,
      #                         y=y, group=grp, color=grp, fill=grp)) +
      #     geom_ribbon(aes(ymin=y_lower,
      #                     ymax=y_upper), alpha=0.3, linetype=0) +
      #     geom_line() +
      #     geom_line(aes(x=x,
      #                   y=y2),
      #               type=2) +
      #     labs(x="distance (m)",y="prediction",title = names(dfvar)[1], color=names(dfvar)[2], fill=names(dfvar)[2])
      # }
    }
    if (!zoi){
      # if (ncol(dfvar) == 1){

      # data for plotting
      df <- data.frame(x = dfvar[,names(dfvar)[1]], y = pred$mid,
                       y_lower = pred$lower, y_upper = pred$higher, y_mean = pred$mean)
      # data for plotting individual models
      if(indiv_pred) df_indiv <- tidyr::gather(data.frame(x = dfvar[,names(dfvar)[1]], pred_indiv),
                                               key = "model", value = "y", -x)

      plt <- ggplot2::ggplot(df)
      if(class(dfvar[,1]) == "factor") {
        if(indiv_pred) {
          plt <- plt + ggplot2::geom_boxplot(ggplot2::aes(x = x,
                                                          y = y),
                                             data = df_indiv)
        }
      } else {
        if(ci) {
          plt <- plt + ggplot2::geom_ribbon(ggplot2::aes(x = x,
                                                         ymin = y_lower,
                                                         ymax = y_upper),
                                            fill = col_ci,
                                            alpha = alpha_ci)
        } else {
          if(indiv_pred) {
            plt <- plt + ggplot2::geom_line(ggplot2::aes(x = x,
                                                         y = y,
                                                         group = model),
                                            data = df_indiv,
                                            color = col_indiv,
                                            linewidth = linewidth_indiv,
                                            alpha = alpha_indiv)
          }
        }
        if(plot_median)
          plt <- plt + ggplot2::geom_line(ggplot2::aes(x=x,
                                                       y = y),
                                          color = col_mean,
                                          linewidth = linewidth_mean)
        if(plot_mean)
          plt <- plt + ggplot2::geom_line(ggplot2::aes(x=x,
                                                       y = y_mean),
                                          color = col_median,
                                          linewidth = linewidth_median)
      }

      plt <- plt + ggplot2::labs(x = names(dfvar), y = y_lab, title = "")

      # }
      # if (ncol(dfvar)==2){
      #   df <- data.frame(x=newdata[,names(dfvar)[1]], grp=as.factor(newdata[,names(dfvar)[2]]), y=pred$mid, y_lower = pred$lower, y_upper=pred$higher, y2=pred$mean)
      #   plt <- ggplot(df, aes(x=x,
      #                         y=y, group=grp, color=grp, fill=grp)) +
      #     geom_ribbon(aes(ymin=y_lower,
      #                     ymax=y_upper), alpha=0.3, linetype=0) +
      #     geom_line() +
      #     geom_line(aes(x=x,
      #                   y=y2),
      #               type=2) +
      #     labs(x=names(dfvar)[1],y="prediction",title = "", color=names(dfvar)[2], fill=names(dfvar)[2])
      # }
    }

    if (logx) { plt <- plt + ggplot2::scale_x_continuous(trans = 'log10') }
    if (!is.null(ylim)) { plt <- plt + ylim }
    return(plt + ggplot2::theme_minimal())

  } else{

    pred <- cbind(dfvar, pred)
    if(indiv_pred) pred <- cbind(pred, pred_indiv)
    return(pred)

  }
}

