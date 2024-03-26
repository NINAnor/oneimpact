#' Adds ZOI radii to formula
#'
#' Adds multiple radii for a zone of influence variable to a formula, to be fitted in
#' a statistical model. If `predictor_table = TRUE`, it also returns a table with information
#' about each predictor (e.g. whether or not it is a ZOI variable, which radius and shape, to
#' what type of infrastructure it corresponds), to be used in the modeling.
#'
#' @param f `[formula]` \cr A formula to be fitted in a statistical model. Here we do not need
#' the actual name of each term, but all ZOI variables whose radius/scale shall be fitted might
#' be added with a pattern, e.g. "road_traffic_zoiXXXX".
#' @param zoi_radius `[numeric,vector]` \cr A vector of radii/scales for the zone of influence,
#' to be added to the formula.
#' @param type `[character=""]` \cr Shape(s) of the ZOI (or vector of shapes if more than one),
#' to be added to the terms name.
#' @param pattern `[character]` \cr Pattern to be replaced in the formula, corresponding to the ZOI
#' radii or ZOI shapes and radii.
#' @param cumulative `[character=""]{"cumulative", "nearest"}` \cr Default is `NULL`. String to be added to the ZOI terms
#' corresponding on whether the variable represents the ZOI of the nearest feature (`cumulative = "nearest"`)
#' or the cumulative ZOI (`cumulative = "cumulative"`). If `""` (default), the type of ZOI is taken from the
#' variable name, already given in the formula.
#' @param separator `[character(1)="_"]` \cr Separator to be used between type and zoi_radius, when `type`
#' is provided. Default is `"_"`.
#' @param predictor_table `[logical(1)=FALSE]` \cr logical. Whether or not a table should be returned
#' with info of all ZOI radii, shape values, and formula terms, together with info from other non-ZOI predictors.
#'
#' @return A list with both the final `formula` with all ZOI radii and shapes and
#' a table with the predictor information (if `predictor_table = TRUE`).
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
#' # predictor_table - adding only the nearest ZOI
#' f <- case_ ~ strata(step_id_) + land_use + sl_*startpt_roads_XXX + sl_*startpt_cabins_XXX
#' add_zoi_formula(f, zoi_radius = c(1000, 2000, 3000), type = c("exp_decay", "threshold"), pattern = "XXX",
#'                 predictor_table = TRUE, cumulative = "nearest")
#'
#' # predictor_table 2 - adding both nearest and cumulative ZOI
#' f <- case_ ~ strata(step_id_) + sl_*startpt_roads_cumulative_XXX + sl_*startpt_cabins_nearest_XXX
#' add_zoi_formula(f, zoi_radius = c(1000, 2000, 3000), type = c("exp_decay", "threshold"), pattern = "XXX",
#'                 predictor_table = TRUE)
#'
#' # predictor_table 3 - adding the nearest/cumulative metric within the type argument
#' f <- case_ ~ strata(step_id_) + sl_*startpt_roads_XXX + sl_*startpt_cabins_XXX
#' add_zoi_formula(f, zoi_radius = c(1000, 2000, 3000), type = c("nearest_exp_decay", "cumulative_exp_decay"), pattern = "XXX",
#'                 predictor_table = TRUE)
#'
#' @export
add_zoi_formula <- function(f, zoi_radius,
                            type = "",
                            pattern = "XXX",
                            cumulative = "",
                            separator = "_",
                            predictor_table = FALSE) {

  # separate case and strata from formula
  f2 <- extract_response_strata(f, other_vars = T)
  # get other variables
  f3 <- strsplit(f2$other_vars, split="+", fixed = T)[[1]] |>
    sapply(trimws, USE.NAMES = FALSE)
  # add zoi_radius
  if(any(type == "")) {
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

  # for variable names
  # f_names <- as.formula(paste0(". ~ - 1 +", f4))
  # all_covs <- all.vars(f_names)[-1]

  # output
  out <- list(formula = f1, predictor_table = predictor_table)
  # fill predictor_table
  if(predictor_table) {
    # get ZOI vals
    lines_zoi <- grep(pattern, grid_zoi$Var3)
    lines_other <- match(unique(grid_zoi$Var3), grid_zoi$Var3)
    lines_other <- lines_other[!(lines_other %in% lines_zoi)]
    grid_zoi$is_zoi <- 1
    grid_zoi[lines_other,c(1,2,4)] <- NA
    grid_zoi <- grid_zoi[sort(c(lines_zoi, lines_other)),]
    grid_zoi$is_zoi <- ifelse(is.na(grid_zoi$is_zoi), 0, grid_zoi$is_zoi)

    # add variable in formula
    #grid_zoi$variable_zoi <- all_covs#f3
    grid_zoi$term_zoi <- f3
    grid_zoi$Var3 <- sapply(grid_zoi$Var3, function(x, y){ gsub(y, "", x)}, y = pattern)
    names(grid_zoi)[c(1,2,3)] <- c("zoi_radius", "shape", "variable")

    # is cumulative?
    if(cumulative != "") {
      grid_zoi$cumulative <- ifelse(grid_zoi$is_zoi, cumulative, NA_character_)
    } else {
      grid_zoi$cumulative <- NA_character_
      grid_zoi$cumulative <- ifelse(grepl("cumulative", grid_zoi$term_zoi), "cumulative", grid_zoi$cumulative)
      grid_zoi$cumulative <- ifelse(grepl("nearest", grid_zoi$term_zoi), "nearest", grid_zoi$cumulative)
    }
    # remove cumulative from shape, if present
    patt_cum <- "cumulative_|cumulative|_cumulative|_cumulative_|nearest|_nearest|nearest_|_nearest_"
    grid_zoi$shape <- gsub(pattern = patt_cum, replacement = "", grid_zoi$shape)

    out$predictor_table <- grid_zoi[, c("is_zoi", "cumulative", "shape", "zoi_radius", "variable", #"variable_zoi",
                                       "term_zoi")]
  }

  return(out)
}
