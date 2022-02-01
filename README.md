# oneimpact

`oneimpact` provides tools for the assessment of cumulative impacts from multiple infrastructure and land use modification in ecological studies.
The tools use R interface but the main calculations might be run in both R and GRASS GIS. The tools available so far are:

#### Influence measures:

- `calc_influence_nearest`: Calculate the influence from the nearest infrastructure, according to multiple possible 
decay functions and zones of influence.
- `calc_influence_cumulative`: Calculate the cumulative influence of multiple infrastructure, according to multiple possible 
decay functions and zones of influence.
- `calc_influence`: Calculate both the influence of the nearest infrastructure and the cumulative influence, at multiple
scales or zones of influence.

#### Auxiliary functions:

- `util_binarize_GRASS`: Binarize continuous or multi-class categorical rasters within GRASS GIS. Binary maps are used 
as input for cumulative influence and kernel density calculation.

#### Spatial filters:

- `create_filter`: Create filters or weight matrices for neighborhood analysis, according to different decay functions
and parameterized using zones of influence.
- `save_rmfilter`: Saves filters/weight matrices outside R for use within GRASS GIS modules.


#### Support for landscape simulation:

- `set_points`: simulate points in a landscape according to different rules and spatial patterns.

## Installation

To install the development version of the `oneimpact` R package, please use:

```
library(devtools)
devtools::install_github("NINAnor/oneimpact", ref = "HEAD")
```

## See also

The `oneimpact` functions are greatly based on neighborhood analyses made through the
[`terra` package](https://rspatial.org/terra/pkg/index.html) in R and on three GRASS GIS modules:
[`r.mfilter`](https://grass.osgeo.org/grass78/manuals/r.mfilter.html), 
[`r.resamp.filter`](https://grass.osgeo.org/grass78/manuals/r.resamp.filter.html), and 
[`r.neighbors`](https://grass.osgeo.org/grass78/manuals/r.neighbors.html). The connection
between R and GRASS GIS is made through the [`rgrass7`](https://github.com/rsbivand/rgrass) R package.

## Meta

  - Please [report any issues or bugs](https://github.com/NINAnor/oneimpact/issues/new/).
  - License: GPL3
  - Get citation information for `oneimpact` in R doing `citation(package = 'oneimpact')`
  - Contributions are mostly welcome!
