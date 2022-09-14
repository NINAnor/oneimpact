#' Sample area: a polygon vector data
#'
#' Dataset containing the limits of an arbitrary study area in Southern Norway, used
#' for illustrative purposes.
#'
#' @name sample_area.gpkg
#'
#' @examples
#' (s <- system.file("vector/sample_area.gpkg", package = "oneimpact"))
#' sf::st_read(s)
#' # or
#' terra::vect(s)
#'
#' @format A geopackage file. Projected CRS: ETRS89 / UTM zone 33N.
NULL

#' Cabins vector data for the sample area
#'
#' Dataset containing the location of tourist private cabins is Southern
#' Norway, within the sample area for the `oneimpact` package.
#' It corresponds to some specific building types (object_type = "Bygning",
#' byggtyp_nbr = c("161", "162", "163")) form the public N50 dataset.
#' The map was cut for the sample area presented in the `oneimpact` package.
#'
#' @name cabins_sample.gpkg
#'
#' @examples
#' (s <- system.file("vector/cabins_sample.gpkg", package = "oneimpact"))
#' sf::st_read(s)
#' # or
#' terra::vect(s)
#'
#' @format A geopackage file. Projected CRS: ETRS89 / UTM zone 33N. The vector
#' file presents the following columns:
#' \itemize{
#'         \item{cat:} {Line number, corresponding to the original dataset}
#'         \item{buildtype:} {Type of building (code) in the original dataset}
#'         \item{city:} {Code of the municipality where the cabin is located}
#'         \item{value:} {Value 1, to be used for rasterization purposes}
#' }
#'
#' @source \url{https://register.geonorge.no/det-offentlige-kartgrunnlaget/n50-kartdata/ea192681-d039-42ec-b1bc-f3ce04c189ac}
NULL

#' Cabin presence raster data for the sample area
#'
#' Raster data indicating pixels with presence of tourist private cabins in Norway.
#' It corresponds to some specific building types (object_type = "Bygning",
#' byggtyp_nbr = c("161", "162", "163")) form the public N50 dataset.
#' The original data consisted of point vector data and were rasterized with 100m
#' resolution, for the purpose of illustration here. The raster was cut for the
#' sample area presented in the `oneimpact` package.
#'
#' @format A Geotiff file. Projected CRS: ETRS89 / UTM zone 33N.
#' \itemize{
#'         \item{1:} {Presence of cabins}
#'         \item{NA:} {No presence of cabins}
#' }
#'
#' @examples
#' (f <- system.file("raster/cabins_sample.tif", package = "oneimpact"))
#' terra::rast(f)
#'
#' @name cabins_sample.tif
#'
#' @source \url{https://register.geonorge.no/det-offentlige-kartgrunnlaget/n50-kartdata/ea192681-d039-42ec-b1bc-f3ce04c189ac}
NULL

#' Cabin count raster data for the sample area
#'
#' Raster data indicating the number of tourist private cabins per pixel in Norway.
#' It corresponds to some specific building types (object_type = "Bygning",
#' byggtyp_nbr = c("161", "162", "163")) form the public N50 dataset.
#' The original data consisted of point vector data and were rasterized with 100m
#' resolution by counting the number of cabins in each pixel. The raster
#' was cut for the study area presented in the `oneimpact` package.
#'
#' @format A Geotiff file. Projected CRS: ETRS89 / UTM zone 33N.
#'
#' @examples
#' (f <- system.file("raster/cabins_sample_count.tif", package = "oneimpact"))
#' terra::rast(f)
#'
#' @name cabins_sample_count.tif
#'
#' @source \url{https://register.geonorge.no/det-offentlige-kartgrunnlaget/n50-kartdata/ea192681-d039-42ec-b1bc-f3ce04c189ac}
NULL

#' Reindeer area: a polygon vector data for the Setesdal Austhei reindeer herding area
#'
#' Dataset containing the limits of the reindeer management area of Setesdal Austhei
#' in Southern Norway, as defined in Panzacchi et al. (2015). Note this area is slightly
#' larger than the boundaries used for reindeer management in Norway.
#'
#' @name reindeer_area.gpkg
#'
#' @examples
#' (s <- system.file("vector/reindeer_area.gpkg", package = "oneimpact"))
#' sf::st_read(s)
#' # or
#' terra::vect(s)
#'
#' @format A geopackage file. Projected CRS: ETRS89 / UTM zone 33N.
#'
#' @source Panzacchi, M., Van Moorter, B., Strand, O., Loe, L. E., & Reimers, E. (2015). Searching
#' for the fundamental niche using individual-based habitat selection modelling across
#' populations. Ecography, 38(7), 659–669. \url{https://doi.org/10.1111/ecog.01075}
NULL

#' Cabins vector data for the reindeer area
#'
#' Dataset containing the location of tourist private cabins is Southern
#' Norway, within the reindeer management area of Setesdal Austhei.
#' It corresponds to some specific building types (object_type = "Bygning",
#' byggtyp_nbr = c("161", "162", "163")) from the public N50 dataset.
#'
#' @name cabins.gpkg
#'
#' @examples
#' (s <- system.file("vector/cabins.gpkg", package = "oneimpact"))
#' sf::st_read(s)
#' # or
#' terra::vect(s)
#'
#' @format A geopackage file. Projected CRS: ETRS89 / UTM zone 33N. The vector
#' file presents the following columns:
#' \itemize{
#'         \item{gid:} {Line number, corresponding to the original dataset}
#'         \item{buildtype:} {Type of building (code) in the original dataset}
#'         \item{city:} {Code of the municipality where the cabin is located}
#'         \item{value:} {Value 1, to be used for rasterization purposes}
#' }
#'
#' @source \url{https://register.geonorge.no/det-offentlige-kartgrunnlaget/n50-kartdata/ea192681-d039-42ec-b1bc-f3ce04c189ac}
NULL

#' Public roads vector data for the reindeer area
#'
#' Dataset containing the location of large, public roads is Southern
#' Norway, within the reindeer management area of Setesdal Austhei.
#' Retrieved from the Norwegian road dataset Elveg 1.0.
#'
#' @name roads_public.gpkg
#'
#' @examples
#' (s <- system.file("vector/roads_public.gpkg", package = "oneimpact"))
#' sf::st_read(s)
#'
#' @format A geopackage file. Projected CRS: ETRS89 / UTM zone 33N. The vector
#' file presents the following columns:
#' \itemize{
#'         \item{id:} {Line number, corresponding to the original dataset}
#'         \item{name:} {Local name of the road}
#'         \item{publ_priv:} {Whether the road is public or private}
#'         \item{traffic_bin:} {Binary classification of the traffic on the road - high or low}
#'         \item{name_area:} {Name of the reindeer management area where the road is located}
#'         \item{traffic_bin:} {Value 1, to be used for rasterization purposes}
#' }
#'
#' @source \url{https://kartkatalog.geonorge.no/metadata/elveg-20/77944f7e-3d75-4f6d-ae04-c528cc72e8f6}
NULL

#' Private roads vector data for the reindeer area
#'
#' Dataset containing the location of small, private roads is Southern
#' Norway, within the reindeer management area of Setesdal Austhei.
#' Retrieved from the Norwegian road dataset Elveg 1.0.
#'
#' @name roads_private.gpkg
#'
#' @examples
#' (s <- system.file("vector/roads_private.gpkg", package = "oneimpact"))
#' sf::st_read(s)
#'
#' @format A geopackage file. Projected CRS: ETRS89 / UTM zone 33N. The vector
#' file presents the following columns:
#' \itemize{
#'         \item{id:} {Line number, corresponding to the original dataset}
#'         \item{name:} {Local name of the road}
#'         \item{publ_priv:} {Whether the road is public or private}
#'         \item{traffic_bin:} {Binary classification of the traffic on the road - high or low}
#'         \item{name_area:} {Name of the reindeer management area where the road is located}
#'         \item{traffic_bin:} {Value 1, to be used for rasterization purposes}
#' }
#'
#' @source \url{https://kartkatalog.geonorge.no/metadata/elveg-20/77944f7e-3d75-4f6d-ae04-c528cc72e8f6}
NULL
