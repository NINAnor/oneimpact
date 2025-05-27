#----------------------------------------------------------------------------------------
# Sample dataset

#' Sample area: a polygon vector data
#'
#' Dataset containing the limits of an arbitrary study area in Southern Norway, used
#' for illustrative purposes.
#'
#' @name sample_area.gpkg
#' @seealso
#' Maps for the sample area: \cr
#' Cabins: [oneimpact::sample_area_cabins.gpkg], [oneimpact::sample_area_cabins.tif],
#' [oneimpact::sample_area_cabins_count.tif] \cr
#' Roads: [oneimpact::sample_area_roads.gpkg], [oneimpact::sample_area_roads.tif]
#'
#' @examples
#' (f <- system.file("vector/sample_area.gpkg", package = "oneimpact"))
#' sf::st_read(f)
#' # or
#' v <- terra::vect(f)
#' plot(v)
#'
#' @format A geopackage file. Projected CRS: ETRS89 / UTM zone 33N.
NULL

#' Cabins vector data for the sample area
#'
#' Dataset containing the location of tourist private cabins in Southern
#' Norway.
#' It corresponds to some specific building types (object_type = "Bygning",
#' byggtyp_nbr = c("161", "162", "163")) form the public N50 dataset.
#' The map was clipped for the sample area presented in the `oneimpact` package.
#'
#' @name sample_area_cabins.gpkg
#' @seealso
#' Maps for the sample area: \cr
#' Limits of sample area: [oneimpact::sample_area.gpkg] \cr
#' Cabins: [oneimpact::sample_area_cabins.tif], [oneimpact::sample_area_cabins_count.tif] \cr
#' Roads: [oneimpact::sample_area_roads.gpkg], [oneimpact::sample_area_roads.tif]
#'
#' @examples
#' (f <- system.file("vector/sample_area_cabins.gpkg", package = "oneimpact"))
#' sf::st_read(f)
#' # or
#' v <- terra::vect(f)
#' plot(v)
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
#' Cabins corresponds to some specific building types (object_type = "Bygning",
#' byggtyp_nbr = c("161", "162", "163")) form the public N50 dataset.
#' The original data consisted of point vector data and were rasterized with 100m
#' resolution, for the purpose of illustration. The raster was clipped for the
#' sample area presented in the `oneimpact` package.
#'
#' @format A Geotiff file. Projected CRS: ETRS89 / UTM zone 33N.
#' \itemize{
#'         \item{1:} {Presence of cabins}
#'         \item{NA:} {No presence of cabins}
#' }
#'
#' @examples
#' (f <- system.file("raster/sample_area_cabins.tif", package = "oneimpact"))
#' r <- terra::rast(f)
#' plot(r)
#'
#' @name sample_area_cabins.tif
#' @seealso
#' Maps for the sample area: \cr
#' Limits of sample area: [oneimpact::sample_area.gpkg] \cr
#' Cabins: [oneimpact::sample_area_cabins.gpkg], [oneimpact::sample_area_cabins_count.tif] \cr
#' Roads: [oneimpact::sample_area_roads.gpkg], [oneimpact::sample_area_roads.tif]
#'
#' @source \url{https://register.geonorge.no/det-offentlige-kartgrunnlaget/n50-kartdata/ea192681-d039-42ec-b1bc-f3ce04c189ac}
NULL

#' Cabin count raster data for the sample area
#'
#' Raster data indicating the number of tourist private cabins per pixel in Norway.
#' Cabins corresponds to some specific building types (object_type = "Bygning",
#' byggtyp_nbr = c("161", "162", "163")) form the public N50 dataset.
#' The original data consisted of point vector data and were rasterized with 100m
#' resolution by counting the number of cabins in each pixel. The raster
#' was clipped for the study area presented in the `oneimpact` package.
#'
#' @format A Geotiff file. Projected CRS: ETRS89 / UTM zone 33N.
#'
#' @examples
#' (f <- system.file("raster/sample_area_cabins_count.tif", package = "oneimpact"))
#' r <- terra::rast(f)
#' plot(r)
#'
#' @name sample_area_cabins_count.tif
#' @seealso
#' Maps for the sample area: \cr
#' Limits of sample area: [oneimpact::sample_area.gpkg] \cr
#' Cabins: [oneimpact::sample_area_cabins.gpkg], [oneimpact::sample_area_cabins.tif] \cr
#' Roads: [oneimpact::sample_area_roads.gpkg], [oneimpact::sample_area_roads.tif]
#'
#' @source \url{https://register.geonorge.no/det-offentlige-kartgrunnlaget/n50-kartdata/ea192681-d039-42ec-b1bc-f3ce04c189ac}
NULL

#' Road vector data for the sample area
#'
#' Dataset containing the location of public roads in Southern
#' Norway.
#' Retrieved from the Norwegian road dataset Elveg 1.0 and
#' clipped for the study area presented in the `oneimpact` package.
#'
#' @name sample_area_roads.gpkg
#' @seealso
#' Maps for the sample area: \cr
#' Limits of sample area: [oneimpact::sample_area.gpkg] \cr
#' Cabins: [oneimpact::sample_area_cabins.gpkg], [oneimpact::sample_area_cabins.tif],
#' [oneimpact::sample_area_cabins_count.tif] \cr
#' Roads: [oneimpact::sample_area_roads.tif]
#'
#' @examples
#' (f <- system.file("vector/sample_area_roads.gpkg", package = "oneimpact"))
#' sf::st_read(f)
#' # or
#' v <- terra::vect(f)
#' plot(v)
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
#' @source \url{https://kartkatalog.geonorge.no/metadata/elveg/ed1e6798-b3cf-48be-aee1-c0d3531da01a}
NULL

#' Road raster data for the sample area
#'
#' Raster data indicating pixels with public roads in Southern
#' Norway. Rasterized from the vector data from the Norwegian road dataset Elveg 1.0
#' with 100 m resolution and
#' clipped for the study area presented in the `oneimpact` package.
#'
#' @format A Geotiff file. Projected CRS: ETRS89 / UTM zone 33N.
#' \itemize{
#'         \item{1:} {Presence of roads}
#'         \item{NA:} {No presence of roads}
#' }
#'
#' @examples
#' (f <- system.file("raster/sample_area_roads.tif", package = "oneimpact"))
#' r <- terra::rast(f)
#' plot(r)
#'
#' @name sample_area_roads.tif
#' @seealso
#' Maps for the sample area: \cr
#' Limits of sample area: [oneimpact::sample_area.gpkg] \cr
#' Cabins: [oneimpact::sample_area_cabins.gpkg], [oneimpact::sample_area_cabins.tif],
#' [oneimpact::sample_area_cabins_count.tif] \cr
#' Roads: [oneimpact::sample_area_roads.gpkg]
#'
#' @source \url{https://kartkatalog.geonorge.no/metadata/elveg/ed1e6798-b3cf-48be-aee1-c0d3531da01a}
NULL

#' Predictor environmental variables for the Hardangervidda wild reindeer area in Norway
#'
#' Raster data with multiple layers representing predictor environmental variables
#' for the Hardangervidda wild reindeer areas in Norway, used for the vignettes and
#' the illustration example in Niebuhr et al. 2023.
#' Comprises multiple layers: land cover, four components of a principal component
#' analysis which represent bio-geo-climatic variation, and 28 layers representing
#' the cumulative zone of influence (ZOI) and the ZOI of the nearest feature for
#' private cottages and public tourist resources, with radii from 100 m to 20 km.
#' The raster has 500 m resolution.
#'
#' @format A Geotiff file. Projected CRS: ETRS89 / UTM zone 33N.
#' \itemize{
#'         \item{norway_pca_klima_axis1-4}{Components 1 to 4 from a principal component analysis
#'         representing bio-geo-climatic variation in Norway, from Bakkestuen et al. 2008.
#'         PCAs 1 to 4 represent, respectively, continentality, altitude, terrain ruggedness, and solar radiation.
#'         More information in Niebuhr et al. 2023. PCAs 1 and 2 also present quadratic layers.}
#'         \item{NORUTreclass} {Land use and land cover classes from NORUT, reclassified as in Niebuhr et al. 2023.}
#'         \item{private_cabins_cumulative_exp_decay_XXX} {Cumulative zone of influence of private cabins at each location,
#'   with exponential decay shape, and radii defined by XXX (from 100 to 20000m).}
#'         \item{private_cabins_nearest_exp_decay_XXX} {Zone of influence of the nearest private cabin at each location,
#'   with exponential decay shape, and radii defined by XXX (from 100 to 20000m).}
#'         \item{public_cabins_high_cumulative_exp_decay_XXX} {Cumulative zone of influence of public resorts at each location,
#'   with exponential decay shape, and radii defined by XXX (from 100 to 20000m).}
#'         \item{public_cabins_high_nearest_exp_decay_XXX} {Zone of influence of the nearest public resort at each location,
#'   with exponential decay shape, and radii defined by XXX (from 100 to 20000m).}
#' }
#'
#' @examples
#' (f <- system.file("raster/rast_predictors_hardanger_500.tif", package = "oneimpact"))
#' r <- terra::rast(f)
#' plot(r)
#'
#' @name rast_predictors_hardanger_500.tif
#' @seealso
#' Data for RSF analysis: [oneimpact::reindeer_rsf]
#'
#' @source Niebuhr, B. B., Van Moorter, B., Stien, A., Tveraa, T., Strand, O., Langeland, K.,
#' Sandström, P., Alam, M., Skarin, A., & Panzacchi, M. (2023). Estimating the cumulative impact
#' and zone of influence of anthropogenic features on biodiversity.
#' Methods in Ecology and Evolution. https://doi.org/10.1111/2041-210X.14133
#'
#' Bakkestuen, V., Erikstad, L., & Halvorsen, R. (2008). Step-less models for regional
#' environmental variation in Norway. Journal of Biogeography, 35(10), 1906–1922.
#' https://doi.org/10.1111/j.1365-2699.2008.01941.x
NULL

#----------------------------------------------------------------------------------------
# Reindeer dataset

#' Reindeer area: a polygon vector data for the Setesdal Austhei reindeer herding area
#'
#' Dataset containing the limits of the reindeer management area of Setesdal Austhei
#' in Southern Norway, as defined in Panzacchi et al. (2015). Note this area is slightly
#' larger than the boundaries used for reindeer management in Norway.
#'
#' @name reindeer_area.gpkg
#'
#' @examples
#' (f <- system.file("vector/reindeer_area.gpkg", package = "oneimpact"))
#' sf::st_read(f)
#' # or
#' v <- terra::vect(f)
#' plot(v)
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
#' @name reindeer_cabins.gpkg
#'
#' @examples
#' (f <- system.file("vector/reindeer_cabins.gpkg", package = "oneimpact"))
#' sf::st_read(f)
#' # or
#' v <- terra::vect(f)
#' plot(v)
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
#' @name reindeer_roads_public.gpkg
#'
#' @examples
#' (f <- system.file("vector/reindeer_roads_public.gpkg", package = "oneimpact"))
#' v <- sf::st_read(f)
#' plot(v[1])
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
#' @name reindeer_roads_private.gpkg
#'
#' @examples
#' (f <- system.file("vector/reindeer_roads_private.gpkg", package = "oneimpact"))
#' v <- sf::st_read(f)
#' plot(v[1])
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
