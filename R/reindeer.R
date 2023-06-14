#' GPS positions for wild reindeer in Norway.
#'
#' @description
#' A dataset containing GPS positions for wild reindeer (*Rangifer tarandus tarandus*) in
#' Setesdal Austhei management area, southern Norway. It includes data from 9 female
#' individuals collected between 2007 and 2010. The coordinates are presented in
#' coordinate reference system ETRS89/UTM 33N (EPSG:25833).
#'
#' Dataset taken from the supplementary data from the paper
#' Cagnacci, F., Focardi, S., Ghisla, A., van Moorter, B., Merrill, E.H., Gurarie,
#' E., Heurich, M., Mysterud, A., Linnell, J., Panzacchi, M., May, R., Nygård, T.,
#' Rolandsen, C. and Hebblewhite, M. (2016), How many routes lead to migration?
#' Comparison of methods to assess and characterize migratory movements. J Anim Ecol,
#' 85: 54-68. https://doi.org/10.1111/1365-2656.12449.
#'
#' @format A data frame with 33,202 rows and 6 variables:
#' \describe{
#'   \item{x}{GPS relocations expressed as UTM easting coordinates, in ETRS89/UTM 33N}
#'   \item{y}{GPS relocations expressed as UTM northing coordinates, in ETRS89/UTM 33N}
#'   \item{t}{Timestamp of GPS relocations, in UTC}
#'   \item{original_animal_id}{Animal name assigned at captures}
#'   \item{animal_year_id}{Unique individual identifier in the population, and sample unit- i.e. data from one animal in one year}
#'   \item{sex}{Sex of the individual animal}
#' }
#'
#' @source Panzacchi, M., Van Moorter, B., Strand, O., Loe, L. E., & Reimers, E. (2015). Searching
#' for the fundamental niche using individual-based habitat selection modelling across
#' populations. Ecography, 38(7), 659–669. \url{https://doi.org/10.1111/ecog.01075}
#'
#' Cagnacci, F., Focardi, S., Ghisla, A., van Moorter, B., Merrill, E.H., Gurarie,
#' E., Heurich, M., Mysterud, A., Linnell, J., Panzacchi, M., May, R., Nygård, T.,
#' Rolandsen, C. and Hebblewhite, M. (2016), How many routes lead to migration?
#' Comparison of methods to assess and characterize migratory movements. J Anim Ecol,
#' 85: 54-68. \url{https://doi.org/10.1111/1365-2656.12449}.
#'
#' Dataset: \url{https://datadryad.org/stash/dataset/doi:10.5061/dryad.rg0v3}
#'
#' @usage data(reindeer)
#'
#' @examples
#' library(tibble)
#' data("reindeer")
#'
"reindeer"

#' Annotated data of wild reindeer in Norway, prepared for step-selection analysis.
#'
#' @description
#' A dataset containing GPS positions for wild reindeer (*Rangifer tarandus tarandus*) in
#' Setesdal Austhei management area, southern Norway, organized in a use-availability design
#' for step-selection analysis (SSA). It includes data from 9 female
#' individuals collected between 2007 and 2010. The coordinates are presented in
#' coordinate reference system ETRS89/UTM 33N (EPSG:25833).
#' The data set was regularized with 3h data and 10 random steps were created for
#' each used step. It was then annotated with rasters on the zone of influence
#' of roads and cabins with radii 1000 - 3000, both at the start and the end
#' point of a step.
#'
#' For more information on the original data, see ?[oneimpact::reindeer].
#'
#' @format A data frame with 31,735 rows and 28 variables:
#' \describe{
#'   \item{x}{GPS relocations expressed as UTM easting coordinates, in ETRS89/UTM 33N}
#'   \item{y}{GPS relocations expressed as UTM northing coordinates, in ETRS89/UTM 33N}
#'   \item{t}{Timestamp of GPS relocations, in UTC}
#'   \item{original_animal_id}{Animal name assigned at captures}
#'   \item{animal_year_id}{Unique individual identifier in the population, and sample unit- i.e. data from one animal in one year}
#'   \item{sex}{Sex of the individual animal}
#'   \item{to be completed}{to be completed}
#' }
#'
#' @source Panzacchi, M., Van Moorter, B., Strand, O., Loe, L. E., & Reimers, E. (2015). Searching
#' for the fundamental niche using individual-based habitat selection modelling across
#' populations. Ecography, 38(7), 659–669. \url{https://doi.org/10.1111/ecog.01075}
#'
#' Cagnacci, F., Focardi, S., Ghisla, A., van Moorter, B., Merrill, E.H., Gurarie,
#' E., Heurich, M., Mysterud, A., Linnell, J., Panzacchi, M., May, R., Nygård, T.,
#' Rolandsen, C. and Hebblewhite, M. (2016), How many routes lead to migration?
#' Comparison of methods to assess and characterize migratory movements. J Anim Ecol,
#' 85: 54-68. \url{https://doi.org/10.1111/1365-2656.12449}.
#'
#' @usage data(reindeer_annotated)
#'
#' @examples
#' library(tibble)
#' data("reindeer_annotated")
#'
"reindeer_annotated"
