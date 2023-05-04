#' Adds ZOI radii to formula
#'
#' Adds multiple radii for zone of influence variable to a formula, to be fitted in
#' a statistical model.
#'
#' @examples
#' f <- case_ ~ strata(step_id_) + sl_*startpt_roadsXXX + sl_*startpt_cabinsXXX
#' add_zoi_formula(f, zoi_radius = c(1000, 2000, 3000), pattern = "XXX")
#'
#' @export
add_zoi_formula <- function(f, zoi_radius, pattern = "XXX") {

  # separate case and strata from formula
  f2 <- extract_case_strata(f, other_vars = T)
  # get other variables
  f3 <- strsplit(f2$other_vars, split="+", fixed = T)[[1]]
  # add zoi_radius
  f3 <- unique(apply(expand.grid(f3, zoi_radius), 1, function(x, y){gsub(y, x[2], x[1])}, y = pattern))
  # re-build formula with all terms
  f3 <- paste0(f3, collapse="+")
  f1 <- paste0(f2$case, " ~ ", "strata(", f2$strata, ") + ", f3)
  f1 <- as.formula(f1)

  return(f1)
}
