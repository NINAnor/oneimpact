#' Plot responses from a bag of models
#'
#' This function takes a bag of models (`x`) and a set of new data (`dfvar`) with variation for one or more
#' specific predictor variables to predict and plot the predictions from the bag. One can either plot only the
#' mean or (weighted) median response for specific preditor variables, and possibly also the
#' confidence interval, computed from the weighted quantiles of the prediction. All other variables
#' are kept constant, as defined by the `baseline` parameter.
#'
#' The function `plot_response` uses the `bag_predict` to produce the predictions.
#'
#' @param x `[bag,list]` \cr A bag of models, resulting from a call to [oneimpact::bag_models()].
#' @param dfvar `[data.frame]` \cr A data.frame with the values of the variables one wants to vary.
#' All other variables are set to their mean or median (this is set by the parameter `baseline`).
#' The column names of the dataframe might correspond exactly to the model covariates or to
#' parts of that (for instance, "roads_paved_" to refer to all ZOI variables related to paved roads).
#' @param data `[data.frame]` \cr The original data used for model fitting. Used only for
#' taking the categories of the categorical variables.
#' @param type `[character(1)="linear"]{"linear", "exponential", "logit", "cloglog"}` \cr Type of response.
#' Might be `"linear"` (default), `"exponential"`, `"logit"`, and `"cloglog"`.
#' @param zoi_shape `[character(1)="linear"]{"exp_decay", "gaussian_decay", "linear_decay", "threshold_decay"}` \cr
#' Shape of the ZOI. Necessary to be specified to represent correctly the estimated ZOI of the
#' preditor variables.
#' @param wq_probs `[vector,numeric(3)=c(0.025, 0.5, 0.975)]` \cr a three element with lower, mid, and higher weighted quantiles to be computed
#' @param ci Should variation or confidence intervals be plotted?
#'
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

  # baselines for plotting and predicting!! mean? median?
  # define the baseline
  bs <- baseline
  if (baseline[1] == "median"){
    baseline <- x$data_summary[rownames(x$data_summary) == "50%",]
  } else {
    if(baseline[1] == "mean") {
      baseline <- x$data_summary[rownames(x$data_summary) == "mean",]
    } else {
      if (baseline[1] == "zero"){
        baseline <- x$data_summary[5,]
        baseline[1, 1+which(x$numeric_covs)] <- 0 # zero for numeric ones
        baseline[1, 1+which(!x$numeric_covs)] <- sapply(unname(which(!x$numeric_covs)),
                                                        function(z) levels(data[,names(x$numeric_covs[z])])[1])
      } else {
        stop(paste0("Invalid 'baseline' parameter: ", baseline, ". Pleas re-define."))
      }
    }
  }

  # is the variable a ZOI variable?
  if (!zoi){

    # new data
    newdata <- baseline
    newdata <- newdata[rep(1, nrow(dfvar)),]
    newdata[,names(dfvar)] <- dfvar

  } else{

    # ZOI radii
    zois <- names(x$data_summary)[grep(names(dfvar)[1], names(x$data_summary))]
    # zois <- as.numeric(gsub(names(dfvar)[1], "", zois))
    zoi_radii <- as.numeric(gsub("\\D", "", zois))

    # compute ZOI for the intended distances
    dfvar2 <- dfvar[,rep(1, length(zoi_radii)), drop = FALSE]
    is_cumulative <- grepl(pattern = which_cumulative, zois)
    dfvar2 <- as.data.frame(do.call("cbind", lapply(c(1:ncol(dfvar2)), function(i) {
      n_feat <- ifelse(is_cumulative[i], n_features, 1)
      n_feat*oneimpact::dist_decay(dfvar2[,i], radius = zoi_radii[i], type = zoi_shape) })))
    names(dfvar2) <- zois#paste0(names(dfvar)[1], zoi_radii)

    # cumulative vars
    #dfvar2[,grepl("cumulative", colnames(dfvar2))]*100

    # new data
    newdata <- baseline
    newdata <- newdata[rep(1, nrow(dfvar)),]
    newdata[,names(dfvar2)] <- dfvar2
    #if (ncol(dfvar)==2){newdata[,names(dfvar)[2]] <- dfvar[,2]}
  }

  # make sure prediction works even if categorical variables are constant
  # get that from the summary instead of dat, where from? maybe a new object
  for(i in which(!x$numeric_covs)) {
    if(class(data[,names(x$numeric_covs[i])]) == "factor") {
      newdata[,i+1] <- factor(newdata[,i+1], levels = levels(data[,names(x$numeric_covs[i])])) # factor
    } else {
      newdata[,i+1] <- factor(newdata[,i+1], levels = sort(unique(data[,names(x$numeric_covs[i])]))) # character
    }
  }

  # predict for new data set

  # predict only for those variables/coefs of interest for ZOI variables
  if(zoi) include <- colnames(dfvar) else include = "all"
  # predict
  pred <- bag_predict(x, newdata, type = type, wMean = T, wq_probs = wq_probs, include = include)

  names(pred) <- c("lower", "mid", "higher", "mean")
  # predict for individual models
  if(indiv_pred) {
    pred_indiv <- bag_predict(x, newdata, type = type, wMean = F, wq_probs = NULL, include = include)[, x$weights > 0]
    # order(x$weights[x$weights > 0])
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
        df <- data.frame(x = newdata[,names(dfvar)[1]], y = pred$mid,
                         y_lower = pred$lower, y_upper = pred$higher, y_mean = pred$mean)
        # data for plotting individual models
        if(indiv_pred) df_indiv <- tidyr::gather(data.frame(x = newdata[,names(dfvar)[1]], pred_indiv),
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

    if (logx) { plt <- plt + scale_x_continuous(trans = 'log10') }
    if (!is.null(ylim)) { plt <- plt + ylim }
    return(plt + ggplot2::theme_minimal())

  } else{

    pred <- cbind(dfvar, pred)
    return(pred)

  }
}

# plot_response_multi_issf_lines <- function(x, dfvar, type=c("linear", "exponential")[2], prop_best=0.5, baseline=NULL, zoi=F, ggplot=T, logx=F, ylim=NULL){
#   if (is.null(baseline)){
#     baseline <- x$data_summary[3,]
#   }
#   if (baseline[1]=="zero"){
#     baseline <- x$data_summary[3,]
#     baseline[1,] <- 0
#   }
#   if (!zoi){
#     newdata <- baseline
#     newdata <- newdata[rep(1, nrow(dfvar)),]
#     newdata[,names(dfvar)] <- dfvar
#   } else{
#     rescaling <- list()
#     rescaling$point <- data.frame(zoi=c(100, 250, 500, 1000, 2500, 5000, 10000),
#                                   scale=c(1.000000, 0.147929, 0.040000, 0.010000, 0.001600, 0.000400, 0.000100))
#
#     rescaling$line <- data.frame(zoi=c(100, 250, 500, 1000, 2500, 5000, 10000),
#                                  scale=c(1.00000000, 0.53846154, 0.29066667, 0.14418182, 0.05788062, 0.02895027, 0.01448610))
#
#     zois <- names(x$data_summary)[grep(names(dfvar)[1], names(x$data_summary))]
#     zois <- as.numeric(gsub(names(dfvar)[1], "", zois))
#     dfvar2 <- dfvar[,rep(1, length(zois))]
#     scal <- rescaling$point$scale[match(zois, rescaling$point$zoi)]
#     dfvar2 <- as.data.frame(do.call("cbind", lapply(c(1:ncol(dfvar2)), function(i,dfvar2,zois,scal){apply(cbind((zois[i]-dfvar2[,i])/zois[i], 0), 1, max)*scal[i]}, dfvar2=dfvar2, zois=zois, scal=scal)))
#     names(dfvar2) <- paste0(names(dfvar)[1], zois)
#
#     newdata <- baseline
#     newdata <- newdata[rep(1, nrow(dfvar)),]
#     newdata[,names(dfvar2)] <- dfvar2
#     if (ncol(dfvar)==2){newdata[,names(dfvar)[2]] <- dfvar[,2]}
#   }
#
#   pred_mean <- predict_multi_issf(x, newdata, type=type, wMean=T)
#   pred <- predict_multi_issf(x, newdata, type=type, wMean=F)
#   incl <- x$weights > quantile(x$weights, props=(1-prop_best/ncol(pred)))
#   pred <- pred[,incl]
#   dm <- dim(pred)
#
#   if (ggplot){
#     require(ggplot2)
#     if (zoi){
#       if (ncol(dfvar)==1){
#         preddf <- data.frame(x=rep(dfvar[,1], dm[2]), y=as.vector(pred), btstrp=rep(c(1:dm[2]), each=dm[1]))
#         df <- data.frame(x=dfvar[,1], y=pred_mean[,1])
#         plt <- ggplot(preddf) +
#           geom_line(aes(x=x, y=y, group=btstrp), color='blue',alpha=5/dm[2]) +
#           geom_line(data=df, aes(x=x, y=y), color='green') + ylim(min(df$y), max(df$y)) +
#           labs(x="distance (m)",y="prediction",title = names(dfvar)[1])
#       }
#       if (ncol(dfvar)==2){
#         df <- data.frame(x=dfvar[,1], grp=as.factor(newdata[,names(dfvar)[2]]), y=pred$mid, y_lower = pred$lower, y_upper=pred$higher, y2=pred$mean)
#         plt <- ggplot(df, aes(x=x,
#                               y=y, group=grp, color=grp, fill=grp)) +
#           geom_ribbon(aes(ymin=y_lower,
#                           ymax=y_upper), alpha=0.3, linetype=0) +
#           geom_line() +
#           geom_line(aes(x=x,
#                         y=y2),
#                     type=2) +
#           labs(x="distance (m)",y="prediction",title = names(dfvar)[1], color=names(dfvar)[2], fill=names(dfvar)[2])
#       }
#     }
#     if (!zoi){
#       if (ncol(dfvar)==1){
#         preddf <- data.frame(x=rep(newdata[,names(dfvar)[1]], dm[2]), y=as.vector(pred), btstrp=rep(c(1:dm[2]), each=dm[1]))
#         df <- data.frame(x=newdata[,names(dfvar)[1]], y=pred_mean[,1])
#         plt <- ggplot(preddf) +
#           geom_line(aes(x=x, y=y, group=btstrp), color='blue',alpha=5/dm[2]) +
#           geom_line(data=df, aes(x=x, y=y), color='green') + ylim(min(df$y), max(df$y)) +
#           labs(x=names(dfvar),y="prediction",title = "")
#       }
#       if (ncol(dfvar)==2){
#         df <- data.frame(x=newdata[,names(dfvar)[1]], grp=as.factor(newdata[,names(dfvar)[2]]), y=pred$mid, y_lower = pred$lower, y_upper=pred$higher, y2=pred$mean)
#         plt <- ggplot(df, aes(x=x,
#                               y=y, group=grp, color=grp, fill=grp)) +
#           geom_ribbon(aes(ymin=y_lower,
#                           ymax=y_upper), alpha=0.3, linetype=0) +
#           geom_line() +
#           geom_line(aes(x=x,
#                         y=y2),
#                     type=2) +
#           labs(x=names(dfvar)[1],y="prediction",title = "", color=names(dfvar)[2], fill=names(dfvar)[2])
#       }
#     }
#     if (logx){plt <- plt + scale_x_continuous(trans='log10')}
#     if (!is.null(ylim)){plt <- plt + ylim}
#     print(plt)
#   } else{
#     pred <- cbind(dfvar, pred)
#     return(pred)
#   }
# }
