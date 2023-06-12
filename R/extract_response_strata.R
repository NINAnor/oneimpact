#' Separates the reponse variable, the strata variable (), and other variables in a statistical formula
#'
#' @param f Formula for model fitting.
#' @param other_vars logical. If true, all the (other) explanatory variables are also
#' appended to the output, separated from the responde and strata variables.
#'
#' @return A `data.frame` with the response variable (`case`), the stratum variable
#' (`strata`), and possibly all explanatory variables, if `other_vars = TRUE`.
#'
#' @examples
#' # formula with strata
#' f <- formula(use ~ exp1 + exp2 + exp1:exp2 + strata(id))
#' extract_response_strata(f)
#' extract_response_strata(f, other_vars = TRUE)
#'
#' # formula without strata
#' f2 <- formula(use ~ exp1 + exp2 + exp1:exp2)
#' extract_response_strata(f2, other_vars = TRUE)
#'
#' @export
extract_response_strata <- function(f, other_vars = F) {

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

    if (other_vars){
      tmp  <- gsub(paste0("strata(", x$strata), 'xx', as.character(f)[3], fixed=T)
      tmp  <- gsub("xx) +", '', tmp, fixed=T) |>
        gsub(pattern = "+ xx)", replacement = '', fixed=T) |>
        gsub(pattern = "xx)", replacement = '', fixed=T)
      x$other_vars <- trimws(tmp)
    }

    return(x)
}
