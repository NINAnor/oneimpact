#' Find layers within GRASS GIS with multiple patterns
#'
#' @param list_patterns `[list,character]` \cr List of strings, each corresponding to a pattern
#' to be searched in the names of layers within a GRASS GIS mapset. The patterns are filtered
#' one after the other, in the order they are listed.
#' @param layers_grass `[vector,character]` \cr Vector of strings with the names of maps being
#' assessed, within a GRASS GIS mapset, such as the ones created through the `g.list` module.
#' If `NULL` (default), the list of maps within the GRASS GIS mapset is assessed within the function.
#' @param type `[character(1)="raster"]` \cr Type of layer to be listed within the GRASS GIS mapset
#' (e.g. "raster", "vector"), if `layers_grass` is `NULL`.
#' @param pattern `[character]` \cr Regular expression used to list maps within the GRASS GIS mapset,
#' if `layers_grass` is `NULL`. Default is "*", so all maps of a given `type` are listed.
#' @param mapset `[character(1)="PERMANENT"]` \cr Name of the mapset from which maps are listed,
#' if `layers_grass` is `NULL`. Default is "PERMANENT".
#'
#' @return One or more strings with the names of the maps within the GRASS GIS mapset.
#'
#' @export
util_find_layer_GRASS <- function(list_patterns,
                             layers_grass = NULL,
                             type = "raster",
                             pattern = "*",
                             mapset = "PERMANENT") {

  # get all raster layers from GRASS, if this is not an input
  if(is.null(layers_grass))
    layers_grass <- rgrass7::execGRASS("g.list", type = type,
                                       pattern = pattern,
                                       mapset = mapset) %>%
      attr(., "resOut")

  # apply all filters according to the patterns
  layers_grass_filt <- layers_grass
  for(i in 1:length(list_patterns)) {
    patt <- list_patterns[[i]]
    layers_grass_filt <- layers_grass_filt %>%
      grep(pattern = patt, value = T)
  }

  # return filtered layer
  layers_grass_filt
}
