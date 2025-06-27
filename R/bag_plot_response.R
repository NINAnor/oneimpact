#' Plot responses from a bag of models
#'
#' This function takes a bag of models (`x`) and a set of new data (`dfvar`) with variation for one or more
#' specific predictor variables to predict and plot the predictions from the bag. One can either plot only the
#' mean or (weighted) median response for specific preditor variables, and possibly also the
#' confidence interval, computed from the weighted quantiles of the prediction. All other variables
#' are kept constant, as defined by the `baseline` parameter.
#'
#' The function `plot_response` uses the `predict` to produce the predictions.
#'
#' @param x `[bag,list]` \cr A bag of models, resulting from a call to [oneimpact::bag_models()].
#' @param dfvar `[data.frame]` \cr A `data.frame` with the values of the variables one wants to vary
#' and predict for.
#' All other variables are set to their mean or median, or to zero (this is set by the parameter `baseline`).
#' The column names of the `data.frame` might correspond exactly to the model covariates or to
#' parts of that (for instance, "roads_paved_" to refer to all ZOI variables related to paved roads).
#' @param data `[data.frame]` \cr The original, complete data used for model fitting. Used only for
#' taking the categories of the categorical variables. Irrelevant if there is no categorical variables.
#' @param type `[character(1)="linear"]{"linear", "exponential", "logit", "cloglog"}` \cr Type of response.
#' Might be `"linear"` (default), `"exponential"`, `"logit"`, and `"cloglog"`.
#' @param zoi_shape `[character(1)="linear"]{"exp_decay", "gaussian_decay", "linear_decay", "threshold_decay"}` \cr
#' Shape of the ZOI. Necessary to be specified to represent correctly the estimated ZOI of the
#' preditor variables.
#' @param wq_probs `[vector,numeric(3)=c(0.025, 0.5, 0.975)]` \cr A three element vector with lower,
#' mid, and higher weighted quantiles to be computed.
#' @param ci Should variation or confidence intervals be plotted?
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
                          zoi_vals = c(100, 250, 500, 1000, 2500, 5000, 10000),
                          type_feature = c("point", "line", "area")[1],
                          type_feature_recompute = FALSE,
                          zoi_limit = 0.05,
                          resolution = 100, # resolution for the raster created in create_line_feature_zoi
                          line_value = 1, # value sey to the linear inftastructure raster created in create_line_feature_zoi
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
                              zoi_vals = c(100, 250, 500, 1000, 2500, 5000, 10000),
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

