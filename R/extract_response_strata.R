#' Separates elements in a statistical formula
#'
#' This function separates the response variable, the strata variable,
#' and the covariates/predictors in a statistical formula.
#'
#' @param f `[formula]` \cr Formula for model fitting.
#' @param covars `[logical(1)=FALSE]` \cr logical. If true, all the (other)
#' covariates/explanatory variables in the formula are also
#' appended to the output, separated from the response and strata variables.
#'
#' @return A `list` with strings representing the response variable (`response`),
#' the stratum variable (`strata`),
#' and possibly all explanatory variables (`covars`), if `covars = TRUE`, as shown
#' in the formula.
#'
#' @examples
#' # formula with strata
#' f <- formula(use ~ exp1 + exp2 + exp1:exp2 + strata(id))
#' extract_response_strata(f)
#' extract_response_strata(f, covars = TRUE)
#'
#' # formula without strata
#' f2 <- formula(use ~ exp1 + exp2 + exp1:exp2)
#' extract_response_strata(f2, covars = TRUE)
#'
#' @export
extract_response_strata <- function(f, covars = FALSE) {

    # create list of variables
    x <- list()
    # transform formula into terms
    tmp <- as.character(f)
    # get the response variable - case 1/0
    x$response  <- tmp[2]

    # get the variable that corresponds to strata
    tmp2 <- grepl("strata", tmp[3])
    if(tmp2 == FALSE) {
      x$strata <- ""
    } else {
      tmp <- substr(tmp[3], gregexpr("strata", tmp[3])[[1]][1]+7, nchar(tmp[3]))
      tmp <- substr(tmp, 1, gregexpr(")", tmp)[[1]][1]-1)
      x$strata <- tmp
    }

    if (covars){
      tmp  <- gsub(paste0("strata(", x$strata), 'xx', as.character(f)[3], fixed=T)
      tmp  <- gsub("xx) +", '', tmp, fixed=T) |>
        gsub(pattern = "+ xx)", replacement = '', fixed=T) |>
        gsub(pattern = "xx)", replacement = '', fixed=T)
      x$covars <- trimws(tmp)
    }

    return(x)
}
