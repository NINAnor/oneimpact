#'---
#' title: oneimpact - Next steps
#' author: "Bernardo Niebuhr"
#' date: "`r format(Sys.time(), '%Y-%m-%d')`"
#' output:
#'   pdf_document: default
#'   html_document: default
#'---
#'
#' 1. `calc_influence_nearest`:
#' - Change `transform` to `type` and include euclidean as the default option
#' - Add support for multiple zoi input values
#' - example in GRASS creating a grass db
#'
#' 2. `calc_influence_cumulative`:
#' - Document parameters and al
#' - implement function within GRASS - r.neighbors and r.resamp.filter
#' - implement CUT extent
#' - make example for GRASS
#' - adapt example creating a grass db
#' - test time for different implementations in GRASS (r.resamp.filter, r.mfilter)
#' - not so important: implement Gaussian with parameters based on zoi
#'
#' 3. `create_filter`:
#' - Document parameters and all
#' - implement Gaussian filter parameterized on zoi
#'
#' 4. All functions:
#' - Possible: standardize parameters (method, type, transform, including "euclidean" as an option for nearest)
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
#'
#' 9. Re-organize examples for GRASS - for reproducibility, one needs to create a grass db, insert the input maps,
#' then run the code.

