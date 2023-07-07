#' Plot responses from a bag of models
#'
#' @param x `[list]` \cr A bag of models, resulting from a call to [oneimpact::bag_models()].
#' @param dfvar `[data.frame]` \cr A data.frame with the values of the variables one wants to vary.
#' All other variables are set to their mean or median (this is set by the parameter `baseline`).
#' The column names of the dataframe might correspond exactly to the model covariates or to
#' parts of that (for instance, "roads_paved_" to refer to all ZOI variables related to paved roads).
#' @param data `[data.frame]` \cr The original data used for model fitting. Used only for
#' taking the categories of the categorical variables.
#' @param type `[character(1)="linear"]{"linear", "exponential"}`
#' @param wQ_probs a three element with lower, mid, and higher weighted quantiles to be computed
#' @param ci Should variation or confidence intervals be plotted?
#'
#' @export
plot_response <- function(x,
                          dfvar,
                          data,
                          type = c("linear", "exponential")[2],
                          zoi_shape = c("exp_decay", "gaussian_decay", "linear_decay", "threshold_decay")[1],
                          which_cumulative = "cumulative",
                          ci = TRUE,
                          wQ_probs = c(0.025, 0.5, 0.975),
                          baseline = c("median", "mean", "zero")[1],
                          zoi = FALSE,
                          zoi_vals = c(100, 250, 500, 1000, 2500, 5000, 10000),
                          ggplot = T,
                          plot_mean = TRUE,
                          plot_median = TRUE,
                          n_features = 1,
                          normalize = c(FALSE, "mean", "median", "ci")[1],
                          logx=F,
                          ylim=NULL){

  UseMethod("plot_response")

}


#' @export
plot_response.bag <- function(x,
                              dfvar,
                              data,
                              type = c("linear", "exponential")[2],
                              zoi_shape = c("exp_decay", "gaussian_decay", "linear_decay", "threshold_decay")[1],
                              which_cumulative = "cumulative",
                              ci = TRUE,
                              wQ_probs = c(0.025, 0.5, 0.975),
                              baseline = c("median", "mean", "zero")[1],
                              zoi = FALSE,
                              zoi_vals = c(100, 250, 500, 1000, 2500, 5000, 10000),
                              ggplot = T,
                              plot_mean = TRUE,
                              plot_median = TRUE,
                              n_features = 1,
                              normalize = c(FALSE, "mean", "median", "ci")[1],
                              logx=F,
                              ylim=NULL){

  # baselines for plotting and predicting!! mean? median?
  # define the baseline
  if (baseline[1] == "median"){
    baseline <- x$data_summary[rownames(x$data_summary) == "50%",]
  } else {
    if(baseline[1] == "mean") {
      baseline <- x$data_summary[rownames(x$data_summary) == "mean",]
    } else {
      if (baseline[1] == "zero"){
        baseline <- x$data_summary[5,]
        baseline[1,] <- 0
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
    dfvar2 <- dfvar[,rep(1, length(zoi_radii))]
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
  for(i in which(!x$numeric_covs)) newdata[,i+1] <- factor(newdata[,i+1], levels = sort(unique(data[,names(x$numeric_covs[i])]))) # maybe different condition and specification with levels if it is a factor

  # predict for new data set
  pred <- bag_predict(x, newdata, type = type, wMean = T, wQ_probs = wQ_probs)
  names(pred) <- c("lower", "mid", "higher", "mean")

  if (ggplot){
    if (zoi){
      if (ncol(dfvar) == 1){

        # data for plotting
        df <- data.frame(x = dfvar[,1], y = pred$mid,
                         y_lower = pred$lower, y_upper = pred$higher, y_mean = pred$mean)

        # normalize y axis?
        if(normalize != FALSE) {
          if(normalize == "mean") range_y <- range(df$y_mean[-1]) else
            if(normalize == "median") range_y <- range(df$y[-1]) else
              range_y <- range(df[-1,2:ncol(df)])
            df[2:ncol(df)] <- do.call("cbind", lapply(2:ncol(df), function(i) (df[,i] - range_y[1])/(diff(range_y))))
        }

        # plot
        plt <- ggplot2::ggplot(df)
        if(ci)
          plt <- plt + ggplot2::geom_ribbon(ggplot2::aes(x = x,
                                                         ymin = y_lower,
                                                         ymax = y_upper),
                                            fill='blue', alpha=0.5)
        if(plot_median)
          plt <- plt + ggplot2::geom_line(ggplot2::aes(x=x,
                                                       y = y),
                                          color='blue')
        if(plot_mean)
          plt <- plt + ggplot2::geom_line(ggplot2::aes(x=x,
                                                       y = y_mean),
                                          color='green')

        plt <- plt + ggplot2::labs(x = ifelse(zoi, "Distance (m)", names(dfvar)), y = "Output", title = "")
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
      if (ncol(dfvar) == 1){
        df <- data.frame(x = newdata[,names(dfvar)[1]], y = pred$mid,
                         y_lower = pred$lower, y_upper = pred$higher, y_mean = pred$mean)
        plt <- ggplot2::ggplot(df)
        if(ci)
          plt <- plt + ggplot2::geom_ribbon(ggplot2::aes(x = x,
                                                         ymin = y_lower,
                                                         ymax = y_upper),
                                            fill='blue', alpha=0.5)
        if(plot_median)
          plt <- plt + ggplot2::geom_line(ggplot2::aes(x=x,
                                                       y = y),
                                          color='blue')
        if(plot_mean)
          plt <- plt + ggplot2::geom_line(ggplot2::aes(x=x,
                                                       y = y_mean),
                                          color='green')
        plt <- plt + ggplot2::labs(x = names(dfvar), y = "Output", title = "")

      }
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
    print(plt)

  } else{
    pred <- cbind(dfvar, pred)
    return(pred)
  }
}


plot_response_multi_issf_lines <- function(x, dfvar, type=c("linear", "exponential")[2], prop_best=0.5, baseline=NULL, zoi=F, ggplot=T, logx=F, ylim=NULL){
  if (is.null(baseline)){
    baseline <- x$data_summary[3,]
  }
  if (baseline[1]=="zero"){
    baseline <- x$data_summary[3,]
    baseline[1,] <- 0
  }
  if (!zoi){
    newdata <- baseline
    newdata <- newdata[rep(1, nrow(dfvar)),]
    newdata[,names(dfvar)] <- dfvar
  } else{
    rescaling <- list()
    rescaling$point <- data.frame(zoi=c(100, 250, 500, 1000, 2500, 5000, 10000),
                                  scale=c(1.000000, 0.147929, 0.040000, 0.010000, 0.001600, 0.000400, 0.000100))

    rescaling$line <- data.frame(zoi=c(100, 250, 500, 1000, 2500, 5000, 10000),
                                 scale=c(1.00000000, 0.53846154, 0.29066667, 0.14418182, 0.05788062, 0.02895027, 0.01448610))

    zois <- names(x$data_summary)[grep(names(dfvar)[1], names(x$data_summary))]
    zois <- as.numeric(gsub(names(dfvar)[1], "", zois))
    dfvar2 <- dfvar[,rep(1, length(zois))]
    scal <- rescaling$point$scale[match(zois, rescaling$point$zoi)]
    dfvar2 <- as.data.frame(do.call("cbind", lapply(c(1:ncol(dfvar2)), function(i,dfvar2,zois,scal){apply(cbind((zois[i]-dfvar2[,i])/zois[i], 0), 1, max)*scal[i]}, dfvar2=dfvar2, zois=zois, scal=scal)))
    names(dfvar2) <- paste0(names(dfvar)[1], zois)

    newdata <- baseline
    newdata <- newdata[rep(1, nrow(dfvar)),]
    newdata[,names(dfvar2)] <- dfvar2
    if (ncol(dfvar)==2){newdata[,names(dfvar)[2]] <- dfvar[,2]}
  }

  pred_mean <- predict_multi_issf(x, newdata, type=type, wMean=T)
  pred <- predict_multi_issf(x, newdata, type=type, wMean=F)
  incl <- x$weights > quantile(x$weights, props=(1-prop_best/ncol(pred)))
  pred <- pred[,incl]
  dm <- dim(pred)

  if (ggplot){
    require(ggplot2)
    if (zoi){
      if (ncol(dfvar)==1){
        preddf <- data.frame(x=rep(dfvar[,1], dm[2]), y=as.vector(pred), btstrp=rep(c(1:dm[2]), each=dm[1]))
        df <- data.frame(x=dfvar[,1], y=pred_mean[,1])
        plt <- ggplot(preddf) +
          geom_line(aes(x=x, y=y, group=btstrp), color='blue',alpha=5/dm[2]) +
          geom_line(data=df, aes(x=x, y=y), color='green') + ylim(min(df$y), max(df$y)) +
          labs(x="distance (m)",y="prediction",title = names(dfvar)[1])
      }
      if (ncol(dfvar)==2){
        df <- data.frame(x=dfvar[,1], grp=as.factor(newdata[,names(dfvar)[2]]), y=pred$mid, y_lower = pred$lower, y_upper=pred$higher, y2=pred$mean)
        plt <- ggplot(df, aes(x=x,
                              y=y, group=grp, color=grp, fill=grp)) +
          geom_ribbon(aes(ymin=y_lower,
                          ymax=y_upper), alpha=0.3, linetype=0) +
          geom_line() +
          geom_line(aes(x=x,
                        y=y2),
                    type=2) +
          labs(x="distance (m)",y="prediction",title = names(dfvar)[1], color=names(dfvar)[2], fill=names(dfvar)[2])
      }
    }
    if (!zoi){
      if (ncol(dfvar)==1){
        preddf <- data.frame(x=rep(newdata[,names(dfvar)[1]], dm[2]), y=as.vector(pred), btstrp=rep(c(1:dm[2]), each=dm[1]))
        df <- data.frame(x=newdata[,names(dfvar)[1]], y=pred_mean[,1])
        plt <- ggplot(preddf) +
          geom_line(aes(x=x, y=y, group=btstrp), color='blue',alpha=5/dm[2]) +
          geom_line(data=df, aes(x=x, y=y), color='green') + ylim(min(df$y), max(df$y)) +
          labs(x=names(dfvar),y="prediction",title = "")
      }
      if (ncol(dfvar)==2){
        df <- data.frame(x=newdata[,names(dfvar)[1]], grp=as.factor(newdata[,names(dfvar)[2]]), y=pred$mid, y_lower = pred$lower, y_upper=pred$higher, y2=pred$mean)
        plt <- ggplot(df, aes(x=x,
                              y=y, group=grp, color=grp, fill=grp)) +
          geom_ribbon(aes(ymin=y_lower,
                          ymax=y_upper), alpha=0.3, linetype=0) +
          geom_line() +
          geom_line(aes(x=x,
                        y=y2),
                    type=2) +
          labs(x=names(dfvar)[1],y="prediction",title = "", color=names(dfvar)[2], fill=names(dfvar)[2])
      }
    }
    if (logx){plt <- plt + scale_x_continuous(trans='log10')}
    if (!is.null(ylim)){plt <- plt + ylim}
    print(plt)
  } else{
    pred <- cbind(dfvar, pred)
    return(pred)
  }
}


#### Find rescaling values from densities:
# zoi <- c(100, 250, 500, 1000, 2500, 5000, 10000)
# rescaling <- data.frame(zoi=zoi, scale=NA)
#
# execGRASS("g.region", parameters = list(vector="wild_reindeer_areas_mid", align="houses_10000@p_prodchange_envpoints"))
#
# areas <- readVECT(vname="wild_reindeer_areas")
# rast <- raster(readRAST("houses_10000@p_prodchange_envpoints", plugin=FALSE))
# plot(rast)
# plot(areas, add=T)
#
# #feat <- SpatialPointsDataFrame(coords=data.frame(x=100000,y=6700000), data=data.frame(id="a"))
# feat <- SpatialLinesDataFrame(SpatialLines(list(Lines(Line(cbind(x=c(25000, 175000),y=c(6650000, 6750000))), ID="a"))), data=data.frame(id="A"), match.ID = F)
# #feat <- SpatialPolygonsDataFrame(SpatialPolygons(list(Polygons(list(Polygon(matrix(c(85000, 6690000, 150000, 6733500, 150000, 6710000, 100000, 6680000), ncol=2, byrow=TRUE))), "a"))), data=data.frame(id="A"), match.ID = F)
#
# writeVECT(feat, vname="tmpXXX", v.in.ogr_flags=c("o","t","overwrite"), driver="SQLite")
# execGRASS("v.to.rast", parameters = list(input="tmpXXX", output="tmp_rast", type="line", use="val"), flags=c("overwrite")) # for points
# execGRASS("v.to.rast", parameters = list(input="tmpXXX", output="tmp_rast", type="line", use="val"), flags=c("d", "overwrite")) # for lines
# execGRASS("r.mapcalc", parameters = list(expression="tmp_rast = if(isnull(tmp_rast), 0, tmp_rast)"), flags="overwrite")
#
# for (i in c(1:nrow(rescaling))){
#   scal=zoi[i]
#   execGRASS("r.resamp.filter", parameters = list(input="tmp_rast", output=paste0("tmp_", scal), filter="bartlett", radius=scal), flags="overwrite")
#   rast <- raster(readRAST(vname=paste0("tmp_", scal, "@u_bram.van.moorter"), plugin=FALSE))
#   plot(rast)
#   rescaling$scale[i] <- max(values(rast))
# }

# sc <- scale(dat$public_cabins_high_cumulative_bartlett_100)
# attr(sc)
# dfvar <- data.frame(public_cabins_high_ = c(1:110)*100)
# data <- dat_sc
# str(data$public_cabins_high_cumulative_bartlett_100)
#

