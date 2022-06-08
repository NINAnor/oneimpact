#'---
#' title: oneimpact - Next steps
#' author: "Bernardo Niebuhr"
#' date: "`r format(Sys.time(), '%Y-%m-%d')`"
#' output:
#'   pdf_document: default
#'   html_document: default
#'---
#'
#' 1. `calc_zoi_nearest`:
#' - Add support for multiple zoi input values for GRASS implementation
#'
#' 2. `calc_influence_cumulative`:
#' - implement function within GRASS - r.neighbors and r.resamp.filter
#' - implement CUT extent
#'
#' 3. `create_filter`:
#' - create filters exactly as r.resamp.filter implementation
#'
#' 5. `calc_influence`:
#' - implement for R and GRASS
#'
#' 6. simulate_dist_dens - with terra - not important
#'
#' 7. Implement calculation using buffers around points instead of the whole map
#'
#' 8. Data:
#' - add a layer to the package, for testing (cabins, describe it)


