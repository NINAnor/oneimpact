#' Computes the conditional Boyce index for model evaluation
#'
#' @param x A data.frame with three columns: x, the predicted values; y, the case variable
#' (use vs. available, 1/0); and strat, the stratum.
#'
#' @export
conditionalBoyce <- function(x,
                             method=c("pearson", "kendall", "spearman")[1],
                             plotit = FALSE,
                             errors = TRUE,
                             warnings = TRUE){
  # x: dataframe with x: predicted value, y:use vs. available, strat: stratum
  # https://cran.r-project.org/web/packages/survival/vignettes/concordance.pdf

  # split by strata
  lx <- split(x, x$strat)
  # length of each strata
  rng <- range(unlist(lapply(lx, nrow)))

  # test for problems in sampling
  test <- as.data.frame(table(unlist(lapply(lx, function(x){ which(x$y[order(x$x)]==1)}))))
  if (errors) {
    if (nrow(test) < 3) stop(paste0("Only ", nrow(test), " bins contain observations, the correlation is unreliable."))
  }
  if (warnings){
    if (rng[1] != rng[2]) warning("Not all strata have the same length.")
    if (nrow(test) < (rng[1]-1)) warning(paste0("Only ", nrow(test), " bins contain observations, the correlation is likely not reliable."))
  }

  # compute correlation
  test2 <- data.frame(bins = c(1:rng[2]), frequency = 0)
  test2$frequency[match(as.character(test[,1]), as.character(test2$bins))] <- test$Freq
  b <- cor(test2$bins, test2$frequency, method = method)
  if (plotit) plot(test2)
  return(b)
}

#' Computes concordance indices for model evaluation
#'
#' These functions compute different concordance indices.
#'
#' @references https://cran.r-project.org/web/packages/survival/vignettes/concordance.pdf
#' @seealso [oneimpact::conditionalBoyce]
#'
#' @name concordance_indices
#' @export
somersD <- function(x,
                    errors = TRUE,
                    warnings = TRUE){
  #x: dataframe with x: predicted value, y:use vs. available, strat: stratum
  #https://cran.r-project.org/web/packages/survival/vignettes/concordance.pdf
  lx <- split(x, x$strat)
  C  <- sum(unlist(lapply(lx, function(x){sum(x$x[x$y==1]>x$x[x$y==0])})))
  D  <- sum(unlist(lapply(lx, function(x){sum(x$x[x$y==1]<x$x[x$y==0])})))
  Tx  <- sum(unlist(lapply(lx, function(x){sum(x$x[x$y==1]==x$x[x$y==0])})))
  d <- (C-D)/(C+D+Tx)
  return(d)
}

#' @rdname concordance_indices
#' @export
AUC <- function(x,
                errors = TRUE,
                warnings = TRUE) {
  d <- somersD(x)
  auc <- (d + 1)/2
  auc
}

# Another possible function for presence/absence responses
# Default cost for binomial outcomes in boot::cv.glm
cost <- function(r, pi = 0, na.rm = TRUE) mean(abs(r-pi) > 0.5, na.rm = na.rm)

