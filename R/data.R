#' Study area: a polygon vector data
#'
#' Dataset containing the limits of an arbitrary study area in Southern Norway, used
#' for illustrative purposes.
#'
#' @name study_area.gpkg
#'
#' @examples
#' (s <- system.file("vector/study_area.gpkg", package = "oneimpact"))
#' sf::st_read(s)
#' # or
#' terra::vect(s)
#'
#' @format A geopackage file. Projected CRS: ETRS89 / UTM zone 33N.
NULL

#' Cabins vector data
#'
#' Dataset containing the location of tourist private cabins is Southern
#' Norway, within the study area for the [oneimpact] package.
#' It corresponds to
#' some specific building types (byggtyp_nbr = c("161", "162", "163",
#' "171", "172")) form the public N50 dataset (describe further here).
#'
#' @name cabins_vect.gpkg
#'
#' @examples
#' (s <- system.file("vector/cabins_vect.gpkg", package = "oneimpact"))
#' sf::st_read(s)
#' # or
#' terra::vect(s)
#'
#' @format A geopackage file. Projected CRS: ETRS89 / UTM zone 33N. The vector
#' file presents the following columns:
#' \itemize{
#'         \item{cat:} {Line number, corresponding to the original dataset}
#'         \item{byggtyp_nbr:} {Type of building (code) in the original dataset}
#'         \item{kommune:} {Code of the municipality where the cabin is located}
#'         \item{value:} {Value 1, to be used for rasterization purposes}
#' }
#'
#' @source \url{to_add_here}
NULL

#' Cabins raster data
#'
#' Raster data indicating pixels with presence of tourist private cabins in Norway. It corresponds to
#' some specific building types (object_type = "Bygning", byggtyp_nbr = c("161", "162", "163",
#' "171", "172")) form the public N50 dataset (describe further here). The original data
#' were point vector data and were rasterized with 100m resolution, for the purpose
#' of illustration here. The raster was cut for the study area presented in the [oneimpact]
#' package.
#'
#' @format A Geotiff file. Projected CRS: ETRS89 / UTM zone 33N.
#' \itemize{
#'         \item{1:} {Presence of cabins}
#'         \item{NA:} {No presence of cabins}
#' }
#'
#' @examples
#' (f <- system.file("raster/cabins.tif", package = "oneimpact"))
#' terra::rast(f)
#'
#' @name cabins.tif
#'
#' @source \url{to_complete_here}
NULL

