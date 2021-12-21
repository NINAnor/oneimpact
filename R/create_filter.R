#' Create filters for raster neighborhood analyses
#' 
#' This function creates matrices of weights following different
#' functions to be used in neighborhood analyses for rasters. It is possible
#' to export these matrices as text files, for use with external softwares
#' such as the r.mfilter module within GRASS GIS.
#' 
#' #### IMPLEMENT BARTLETT
#' #### POLSSIBLY: IMPLEMENT IN THE SAME WAY AS FOCAL, WITH INPUT RASTER AS ARGUMENT, POSSIBLY
#' 
#' @inheritParams set_filt_exp_decay
#' 
#' @return A matrix with the weight values.
#' 
#' @example examples/create_filter_example.R
#' 
#' @export
create_filter <- function(zoi = NULL,
                          res = 100,
                          method = c("exp_decay", "step", "bartlett")[1],
                          half_life = NULL,
                          zoi_hl_ratio = 4,
                          min_intensity = 0.01,
                          max_dist = 5000,
                          normalize = FALSE,
                          round_vals = 3,
                          save_txt = FALSE,
                          save_format = "GRASS",
                          save_folder = ".",
                          output = c("CumInf", "Densiy")[1],
                          paralell = TRUE) {
  
  # apply function
  if(method == "exp_decay") {
    parms <- set_filt_exp_decay(zoi = zoi, half_life = half_life,
                                res = res, zoi_hl_ratio = zoi_hl_ratio,
                                min_intensity = 0.01, max_dist = 50000)
  } else {
    if(method == "step") {
      parms <- set_filt_step(zoi = zoi, res = res)
    } else {
      if(method == "bartlett") {
        parms <- set_filt_bartlett(zoi = zoi, res = res)
      }
    }
  }
  
  # get parameters
  zoi <- parms$zoi
  radius_pix <- parms$radius_pix
  size_pix <- parms$size_pix
  
  # create distance matrix
  # distance in pixels to the central cell of the matrix
  dist_mat <- sqrt((matrix(c(1:size_pix), nrow = size_pix, ncol = size_pix, byrow = F) - (radius_pix + 1))^2+
                     (matrix(c(1:size_pix), nrow = size_pix, ncol = size_pix, byrow = T) - (radius_pix + 1))^2)
  # image(dist_mat)
  # plot(terra::rast(dist_mat))
  
  # apply function
  if(method == "exp_decay") {
    dist_mat <- exp(-parms$lambda * dist_mat)
  } else {
    if(method == "step") {
      dist_mat <- 1 * (dist_mat <= radius_pix)
    } else {
      if(method == "bartlett") {
        dist_mat <- pmax((1 + parms$beta * dist_mat), 0)
      }
    }
  }
  # image(dist_mat)
  # plot(terra::rast(dist_mat))
  
  # normalize
  if(normalize)
    dist_mat <- dist_mat/sum(dist_mat[1+radius_pix,])
  # dist_mat2 <- dist_mat/sum(dist_mat)
  #### OR NORMALIZE COMPLETELY??
  # image(dist_mat)
  # plot(terra::rast(dist_mat))
  
  # round decimals
  if(round_vals >= 0) dist_mat <- round(dist_mat, round_vals)
  # image(dist_mat)
  # plot(terra::rast(dist_mat))
  
  if(save_txt) {
    if(output == "CumInf") DIV <- 0 else DIV <- size_pix**2
    if(paralell) TYP <- "P" else TYP <- "S"
    if(save_format == "GRASS") {
       con <- file(paste0(save_folder, "/filter_exp_decay_", zoi, "m.txt"), "w")
       writeLines(paste0('TITLE "filter exp decay ', half_life_input, 'm"'), con = con, sep = "\n", useBytes = FALSE)
       writeLines(paste0("MATRIX ", size_pix), con = con, sep = "\n", useBytes = FALSE)
       for (i in c(1:nrow(dmat))) {
         writeLines(paste0(dmat[i,], collapse=" "), con = con, sep = "\n", useBytes = FALSE)
       }
       writeLines(paste0("DIVISOR ", DIV), con = con, sep = "\n", useBytes = FALSE)
       writeLines(paste0("TYPE ", TYP), con = con, sep = "\n", useBytes = FALSE)
       close(con)
    }
  }
  
  dist_mat
}

#' @param zoi `[numeric(1)]` \cr Zone of Influence (ZoI), in meters. The ZoI is distance or scale up to 
#' which we consider there is effect of an infrastructure or variable for setting the filter.
#' @param half_life `[numeric(1)]` \cr Half life parameter of the exponential decay function, in meters. If NULL,
#' the half life is define in terms of the ZoI and the `zoi_hl_ratio` parameter, which defines the ratio 
#' between the ZoI and the half life. By default, we set this ratio as `zoi/half_life = 4`.
#' The exponent of the exponential decay distance function is defined as `lambda = log(2)/half_life`.
#' @param res `[numeric(1)=100]` \cr Resolution (pixel size) of the filter.
#' @param zoi_hl_ratio `[numeric(1)=6]` \cr Ratio between the ZoI and the half life of the exponential decay
#' distance function. It is used to define the ZoI for the exponential decay function. For instance, if 
#' `half_life = 1000` and `zoi_hl_ratio = 4`, the ZoI will be 4000 m (when the exponential decay decrease to 
#' `0.5**4 = 0.0625`.
#' @param min_intensity `[numeric(1)=0.01]` \cr Minimum intensity of the exponential decay function to
#' define the size (radius) of the window that define the filter.
#' @param max_dist `[numeric(1)=50000]` \cr Maximum size (in meters) to
#' define the size (radius) of the window that define the filter.
set_filt_exp_decay <- function(zoi = NULL,
                               half_life = NULL,
                               res = 100,
                               zoi_hl_ratio = 4,
                               min_intensity = 0.01,
                               max_dist = 50000){
  
  # define zoi or half life, depending on which is given as input
  if(is.null(half_life)) 
    half_life <- zoi/zoi_hl_ratio
  else
    zoi <- half_life * zoi_hl_ratio
  
  # define half life in terms on number of pixels
  half_life <- half_life/res
  # define lambda
  lambda <- log(2)/half_life
  # tmp <- exp(-lambda * c(0:round(half_life*6))/half_life)
  # define radius and size (diameter)
  tmp <- exp(-lambda * c(0:round(2*zoi)))
  radius_pix <- min(which(tmp < min_intensity)[1], round(max_dist/res))
  size_pix <- 2*radius_pix + 1
  
  return(list(zoi = zoi, radius_pix = radius_pix, size_pix = size_pix, lambda = lambda))
}

set_filt_step <- function(zoi, res){
  
  # define radius and size (diameter)
  radius_pix <- floor(zoi/res)
  size_pix <- 2*radius_pix + 1
  
  return(list(zoi = zoi, radius_pix = radius_pix, size_pix = size_pix, lambda = NULL))
}

set_filt_bartlett <- function(zoi, res){

  # define radius and size (diameter)
  radius_pix <- floor(zoi/res)
  size_pix <- 2*radius_pix + 1
  # define beta (beta = -b/a or beta = -1/zoi)
  beta <- -1/radius_pix

  return(list(zoi = zoi, radius_pix = radius_pix, size_pix = size_pix, beta = beta))
}

# ed(100)
# ed(250)
# ed(500)
# ed(1000)
# ed(2500)
# ed(5000, maxdist=10000)
# 
# meanlife=5000
# plot(c(0:20000), exp(-c(0:20000)/meanlife))
# abline(v=meanlife)
# abline(h=0.367879441)

#' 
#' to use with GRASS' r.mfilter
#' 
#' 
#' 

# meanlife=1000
# ed <- function(meanlife, maxdist=5000, res=100){
#   meanlife <- meanlife/res
#   tmp <- exp(-c(0:round(meanlife*5))/meanlife)
#   radius_pix <- min(which(tmp<0.05)[1], round(maxdist/res))
#   dmat <- sqrt((matrix(c(1:(1+2*radius_pix)), nrow=1+2*radius_pix, ncol=1+2*radius_pix, byrow=F)-(radius_pix+1))^2+
#                  (matrix(c(1:(1+2*radius_pix)), nrow=1+2*radius_pix, ncol=1+2*radius_pix, byrow=T)-(radius_pix+1))^2)
#   image(dmat)
#   dmat <- exp(-dmat/meanlife)
#   image(dmat)
#   dmat <- dmat/sum(dmat[1+radius_pix,])
#   image(dmat)
#   dmat <- round(dmat, 3)
#   image(dmat)
#   
#   con <- file(paste0("meandist", meanlife*res, ".txt"), "w")
#   writeLines(paste0("TITLE meandist",meanlife), con = con, sep = "\n", useBytes = FALSE)
#   writeLines(paste0("MATRIX ",(1+2*radius_pix)), con = con, sep = "\n", useBytes = FALSE)
#   for (i in c(1:nrow(dmat))){
#     writeLines(paste0(dmat[i,], collapse=" "), con = con, sep = "\n", useBytes = FALSE)
#   }
#   writeLines(paste0("DIVISOR ",1), con = con, sep = "\n", useBytes = FALSE)
#   writeLines("TYPE P", con = con, sep = "\n", useBytes = FALSE)
#   close(con)
# }
# 
# ed(100)
# ed(250)
# ed(500)
# ed(1000)
# ed(2500)
# ed(5000, maxdist=10000)
# 
# meanlife=5000
# plot(c(0:20000), exp(-c(0:20000)/meanlife))
# abline(v=meanlife)
# abline(h=0.367879441)
