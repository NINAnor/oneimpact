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
#' @usage data(reindeer_ssf)
#'
#' @examples
#' library(tibble)
#' data("reindeer_ssf")
#'
"reindeer_ssf"

#' Annotated data of wild reindeer in Norway, prepared for point resource-selection analysis.
#'
#' @description
#' An data set where GPS positions of wild reindeer (*Rangifer tarandus tarandus*)
#' Hardangervidda management area, southern Norway, were organized in a use-availability design
#' for point resource-selection analysis (RSA) and annotated with environmental data.
#' It includes data from 115 female individuals collected between 2001 and 2019.
#' The data set was regularized with 3h data and 9 random locations were created for
#' each used location. It was then annotated with rasters on the zone of influence
#' of private cabins and public resorts with exponential decay shape and radii
#' varying from 100m to 10km, as well as land cover
#' and bio-geo-climatic variables. Geographical coordinates of the used and random positions
#' were omitted after data annotation.
#'
#' This is part of the data set used for analysis in Niebuhr et al. (2023). ZOI variables
#' with shapes other than the exponential decay were omitted.
#'
#' @format A data frame with 31,735 rows and 28 variables:
#' \describe{
#'   \item{use}{Case in the use-availability setup; 1 represents a used location, 0 represents a random locations.}
#'   \item{private_cabins_cumulative_exp_decay_XXX}{Cumulative zone of influence of private cabins at each location,
#'   with exponential decay shape, and radii defined by XXX (from 100 to 10000m).}
#'   \item{private_cabins_nearest_exp_decay_XXX}{Zone of influence of the nearest private cabin at each location,
#'   with exponential decay shape, and radii defined by XXX (from 100 to 10000m).}
#'   \item{public_cabins_cumulative_exp_decay_XXX}{Cumulative zone of influence of public resorts at each location,
#'   with exponential decay shape, and radii defined by XXX (from 100 to 10000m).}
#'   \item{public_cabins_nearest_exp_decay_XXX}{Zone of influence of the nearest public resort at each location,
#'   with exponential decay shape, and radii defined by XXX (from 100 to 10000m).}
#'   \item{NORUTreclass}{Land use and land cover classes from NORUT, reclassified as in Niebuhr et al. 2023.}
#'   \item{norway_pca_klima_axis1-4}{PCAs 1 to 4 representing bio-geo-climatic variation in Norway,
#'   from Bakkestuen et al. 2008. More information in Niebuht et al. 2023.}
#' }
#'
#' @source Niebuhr, B. B., Van Moorter, B., Stien, A., Tveraa, T., Strand, O., Langeland, K.,
#' Sandström, P., Alam, M., Skarin, A., & Panzacchi, M. (2023). Estimating the cumulative impact
#' and zone of influence of anthropogenic features on biodiversity.
#' Methods in Ecology and Evolution. https://doi.org/10.1111/2041-210X.14133
#'
#' Bakkestuen, V., Erikstad, L., & Halvorsen, R. (2008). Step-less models for regional
#' environmental variation in Norway. Journal of Biogeography, 35(10), 1906–1922.
#' https://doi.org/10.1111/j.1365-2699.2008.01941.x
#'
#' @usage data(reindeer_rsf)
#'
#' @examples
#' library(tibble)
#' data("reindeer_rsf")
#'
"reindeer_rsf"
