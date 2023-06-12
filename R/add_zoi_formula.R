#' Adds ZOI radii to formula
#'
#' Adds multiple radii for a zone of influence variable to a formula, to be fitted in
#' a statistical model.
#'
#' @param f A formula to be fitted in a statistical model.
#' @param zoi_radius `[numeric]` \cr A vector of radii for the zone of influence, to be added to
#' the formula.
#' @param type Shape(s) of the ZOI
#' @param pattern Pattern to be replaced in the formula, corresponding to the ZOI radii or
#' ZOI shapes and radii.
#' @param separator Separator to be used between type and zoi_radius, when `type`
#' is provided. Default is "_".
#' @param grid logical. Whether or not the grid with all ZOI radii and shape values
#' should be returned.
#'
#' @examples
#' # multiple radii
#' f <- case_ ~ strata(step_id_) + sl_*startpt_roadsXXX + sl_*startpt_cabinsXXX
#' add_zoi_formula(f, zoi_radius = c(1000, 2000, 3000), pattern = "XXX")
#'
#' # multiple radii and shapes
#' f <- case_ ~ strata(step_id_) + sl_*startpt_roadsXXX + sl_*startpt_cabinsXXX
#' add_zoi_formula(f, zoi_radius = c(1000, 2000, 3000), type = "_exp_decay", pattern = "XXX")
#' add_zoi_formula(f, zoi_radius = c(1000, 2000, 3000), type = c("_exp_decay", "_threshold"), pattern = "XXX")
#'
#' @export
add_zoi_formula <- function(f, zoi_radius, type = NULL, pattern = "XXX",
                            separator = "_", grid = FALSE) {

  # separate case and strata from formula
  f2 <- extract_response_strata(f, other_vars = T)
  # get other variables
  f3 <- strsplit(f2$other_vars, split="+", fixed = T)[[1]] |>
    sapply(trimws, USE.NAMES = FALSE)
  # add zoi_radius
  if(is.null(type)) {
    grid_zoi <- expand.grid(zoi_radius, type = NA, f3)
    f3 <- unique(apply(grid_zoi, 1,
                       function(x, y){ gsub(y, as.numeric(x[1]), x[3])}, y = pattern))
  } else {
    grid_zoi <- expand.grid(zoi_radius, type, f3)
    f3 <- unique(apply(grid_zoi, 1,
                       function(x, y){ gsub(y, paste0(x[2], separator, as.numeric(x[1])),
                                            x[3])}, y = pattern))
  }
  # re-build formula with all terms
  f4 <- paste0(f3, collapse = " + ")
  # add strata only if strata exists
  f1 <- paste0(f2$response, " ~ ", ifelse(f2$strata == "", "", paste0("strata(", f2$strata, ") + ")), f4)
  # transform into formula
  f1 <- as.formula(f1)

  # output
  out <- list(formula = f1, grid = NA)
  # fill grid
  if(grid) {
    lines_zoi <- grep(pattern, grid_zoi$Var3)
    lines_other <- match(unique(grid_zoi$Var3), grid_zoi$Var3)
    lines_other <- lines_other[!(lines_other %in% lines_zoi)]
    grid_zoi[lines_other,c(1,2)] <- NA
    grid_zoi <- grid_zoi[c(lines_zoi, lines_other),][-3]
    grid_zoi$variable <- f3
    names(grid_zoi)[c(1,2)] <- c("zoi_radius", "shape")
    out$grid <- grid_zoi
  }

  return(out)
}
