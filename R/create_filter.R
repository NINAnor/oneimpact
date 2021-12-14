#' Create filters for raster neighborhood analyses
#' 
#' This function creates matrices of weights following different
#' function to be used in neighborhood analyses for rasters. It is possible
#' to export these matrices as text files, for use with external softwares
#' such as the r.mfilter module within GRASS GIS.
#' 
#' @return A matrix with the weight values.
#' 
#' @example examples/create_filter_example.R
#' 
#' @export

ed <- function(exp_hl = 1000, method = c("exp_decay"), 
               max_dist = 5000, 
               res = 100,
               save_txt = FALSE) {
  
  if(method == "exp_decay") {
    exp_hl <- exp_hl/res
    tmp <- exp(-c(0:round(meanlife*5))/meanlife)
    sz <- min(which(tmp<0.05)[1], round(maxdist/res))
  }
  
  dmat <- sqrt((matrix(c(1:(1+2*sz)), nrow=1+2*sz, ncol=1+2*sz, byrow=F)-(sz+1))^2+
                 (matrix(c(1:(1+2*sz)), nrow=1+2*sz, ncol=1+2*sz, byrow=T)-(sz+1))^2)
  image(dmat)
  dmat <- exp(-dmat/meanlife)
  image(dmat)
  dmat <- dmat/sum(dmat[1+sz,])
  image(dmat)
  dmat <- round(dmat, 3)
  image(dmat)
  
  if(save_txt) {
    con <- file(paste0("meandist", meanlife*res, ".txt"), "w")
    writeLines(paste0("TITLE meandist",meanlife), con = con, sep = "\n", useBytes = FALSE)
    writeLines(paste0("MATRIX ",(1+2*sz)), con = con, sep = "\n", useBytes = FALSE)
    for (i in c(1:nrow(dmat))){
      writeLines(paste0(dmat[i,], collapse=" "), con = con, sep = "\n", useBytes = FALSE)
    }
    writeLines(paste0("DIVISOR ",1), con = con, sep = "\n", useBytes = FALSE)
    writeLines("TYPE P", con = con, sep = "\n", useBytes = FALSE)
    close(con)
  }
  
}

ed(100)
ed(250)
ed(500)
ed(1000)
ed(2500)
ed(5000, maxdist=10000)

meanlife=5000
plot(c(0:20000), exp(-c(0:20000)/meanlife))
abline(v=meanlife)
abline(h=0.367879441)

#' 
#' to use with GRASS' r.mfilter
#' 
#' 
#' 

meanlife=1000
ed <- function(meanlife, maxdist=5000, res=100){
  meanlife <- meanlife/res
  tmp <- exp(-c(0:round(meanlife*5))/meanlife)
  sz <- min(which(tmp<0.05)[1], round(maxdist/res))
  dmat <- sqrt((matrix(c(1:(1+2*sz)), nrow=1+2*sz, ncol=1+2*sz, byrow=F)-(sz+1))^2+
                 (matrix(c(1:(1+2*sz)), nrow=1+2*sz, ncol=1+2*sz, byrow=T)-(sz+1))^2)
  image(dmat)
  dmat <- exp(-dmat/meanlife)
  image(dmat)
  dmat <- dmat/sum(dmat[1+sz,])
  image(dmat)
  dmat <- round(dmat, 3)
  image(dmat)
  
  con <- file(paste0("meandist", meanlife*res, ".txt"), "w")
  writeLines(paste0("TITLE meandist",meanlife), con = con, sep = "\n", useBytes = FALSE)
  writeLines(paste0("MATRIX ",(1+2*sz)), con = con, sep = "\n", useBytes = FALSE)
  for (i in c(1:nrow(dmat))){
    writeLines(paste0(dmat[i,], collapse=" "), con = con, sep = "\n", useBytes = FALSE)
  }
  writeLines(paste0("DIVISOR ",1), con = con, sep = "\n", useBytes = FALSE)
  writeLines("TYPE P", con = con, sep = "\n", useBytes = FALSE)
  close(con)
}

ed(100)
ed(250)
ed(500)
ed(1000)
ed(2500)
ed(5000, maxdist=10000)

meanlife=5000
plot(c(0:20000), exp(-c(0:20000)/meanlife))
abline(v=meanlife)
abline(h=0.367879441)
