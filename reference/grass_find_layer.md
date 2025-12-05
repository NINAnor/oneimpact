# Find layers within GRASS GIS with multiple patterns

Find layers within GRASS GIS with multiple patterns

## Usage

``` r
grass_find_layer(
  list_patterns,
  layers_grass = NULL,
  type = "raster",
  pattern = "*",
  mapset = "PERMANENT"
)
```

## Arguments

- list_patterns:

  `[list,character]`  
  List of strings, each corresponding to a pattern to be searched in the
  names of layers within a GRASS GIS mapset. The patterns are filtered
  one after the other, in the order they are listed.

- layers_grass:

  `[vector,character]`  
  Vector of strings with the names of maps being assessed, within a
  GRASS GIS mapset, such as the ones created through the `g.list`
  module. If `NULL` (default), the list of maps within the GRASS GIS
  mapset is assessed within the function.

- type:

  `[character(1)="raster"]`  
  Type of layer to be listed within the GRASS GIS mapset (e.g. "raster",
  "vector"), if `layers_grass` is `NULL`.

- pattern:

  `[character]`  
  Regular expression used to list maps within the GRASS GIS mapset, if
  `layers_grass` is `NULL`. Default is "\*", so all maps of a given
  `type` are listed.

- mapset:

  `[character(1)="PERMANENT"]`  
  Name of the mapset from which maps are listed, if `layers_grass` is
  `NULL`. Default is "PERMANENT".

## Value

One or more strings with the names of the maps within the GRASS GIS
mapset.
