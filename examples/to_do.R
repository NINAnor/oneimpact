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
#' 
#' - Add support for multiple zoi input values
#' 
#' 2. `calc_influence_cumulative`:
#' - Document parameters and al
#' - implement exp_decay and Bartlett without the need to use create_filter (embed it)
#' - implement function within GRASS
#' - make example for GRASS
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

