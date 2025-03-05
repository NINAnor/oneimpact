#' Create ZOI layer names as strings for data annotation
#'
#' This function uses a vector of variable names and a vector of ZOI radii to create
#' a `data.frame` with all combinations of variable names and zoi radii, to ease
#' listing or accessing them from GRASS GIS or for setting up a model formula.
#'
#' @examples
#' # attaching zoi_radius to the end
#' name <- c("houses", "roads", "railways")
#' radii <- seq(100, 500, by = 100)
#' add_zoi_layers(name, radii)
#'
#' # replacing a pattern by zoi_radius
#' themes <- c("houses_XXX", "private_cabins_XXX", "roads_XXX", "powerlines_XXX", "railways_XXX")
#' maps <- c("houses_XXX@my_mapset", "private_cabins_XXX@my_mapset", "roads_summer_XXX@my_mapset",
#'          "powerlines_XXX@my_mapset", "railway_XXX@my_mapset")
#'
#' add_zoi_layers(layers = maps, zoi_radius = c(100, 500, 1000, 5000, 10000),
#'                name = themes, pattern = "XXX")
#'
#' @export
add_zoi_layers <- function(layers, zoi_radius, name = NULL, pattern = NULL) {

  # options for patterns

  # attach to the end
  if(is.null(pattern)) {
    f <- function(x, y) paste0(x[1], as.numeric(x[2]))
  } else {
    # replace the pattern
    f <- function(x, y) gsub(y, as.numeric(x[2]), x[1])
  }

  # replace layers
  gr <- expand.grid(layers, zoi_radius)
  layers_zoi <- unique(apply(gr, 1, f, y = pattern) )
  zoi_vals <- gr[,2]

  # replace names
  if(!is.null(name)) {
    names_zoi <- unique(apply(expand.grid(name, zoi_radius), 1, f, y = pattern))
    layers_out <- data.frame(name = names_zoi, zoi_radius = zoi_vals, layer = layers_zoi)
  } else {
    layers_out <- data.frame(zoi_radius = zoi_vals, layer = layers_zoi)
  }

  return(layers_out)
}
