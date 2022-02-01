#' Cabins raster data
#'
#' Raster data with locations of tourist private cabins in Norway. It corresponds to
#' some specific building types (object_type = "Bygning", byggtyp_nbr = c("161", "162", "163",
#' "171", "172")) form the public N50 dataset (describe further here). The original data
#' were point vector data and were rasterized with 100m resolution for the purpose
#' of illustration here. Value 1 corresponds to pixels with cabins, value NA corresponds to
#' pixels without cabins. The raster was cut for the study area presented here in the
#' package.
#'
#' @format A RasterLayer object. Projected CRS: ETRS89 / UTM zone 33N.
#'
#' @source \url{}
"cabins"

#' Study area: a polygon vector data
#'
#' Dataset containing the limits of an arbitrary study area in Southern Norway, used
#' for illustrative purposes.
#'
#' @format A sf object. Projected CRS: ETRS89 / UTM zone 33N.
"study area"
