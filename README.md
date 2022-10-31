# oneimpact <img src="man/figures/oneimpact_hex_logo.png" align="right" alt="" width="150" />

<!-- badges: start -->
  [![DOI](https://zenodo.org/badge/453101311.svg)](https://zenodo.org/badge/latestdoi/453101311)
<!--  [![R-CMD-check](https://github.com/NINAnor/oneimpact/workflows/R-CMD-check/badge.svg)](https://github.com/NINAnor/oneimpact/actions) -->
<!-- badges: end -->

`oneimpact` provides tools for the assessment of cumulative impacts of multiple infrastructure and land use modifications in ecological studies.
The tools use R interface but the main calculations might be run in both R and GRASS GIS. The tools available so far are:

### Zone of influence (ZoI) decay functions

- zoi_functions: a set of decay zone of influence functions to characterize different shapes of the ZoI around infrastructure, 
parameterized based on the zone of influence radius. The functions implemented so far are: threshold (`threshold_decay` or `step_decay`),
linear decay (`linear_decay` or `bartlett_decay` or `tent_decay`), exponential decay (`exp_decay`), or Gaussian decay 
(`gaussian_decay` or `half_norm_decay`).
- `plot_zoi1d`: plot ZoI in 1 dimensional space for multiple points infrastructure, using both the ZoI of the nearest
feature and the cumulative ZoI metric.


### Compute zones of influence (ZoI):

- `calc_zoi_nearest`: Calculate the zone of influence from the nearest infrastructure, according to multiple possible 
decay functions and zones of influence radii.
- `calc_zoi_cumulative`: Calculate the cumulative zone of influence of multiple features, according to multiple possible 
decay functions and zones of influence radii.
- `calc_zoi`: Calculate both the the ZoI of the nearest infrastructure and the cumulative ZoI, at multiple
scales or zones of influence radii.

### Spatial filters:

- `create_filter`: Create filters or weight matrices for neighborhood analysis, according to different decay functions
and parameterized using the zone of influence radius.
- `save_filter`: Saves filters/weight matrices outside R for use within GRASS GIS modules.

### Ancillary functions:

- `grass_binarize`: Binarize continuous or multi-class categorical rasters within GRASS GIS. Binary maps may be used 
as input for cumulative zone of influence and kernel density calculation.
- `grass_v2rast_count`: Rasterize a vector files counting the number of features within each pixel of the output
raster. Count rasters may be used as input for cumulative zone of influence and kernel density calculation.

### Support for landscape simulation:

- `set_points`: simulate points in a landscape according to different rules and spatial patterns.

## Installation

To install the development version of the `oneimpact` R package, please use:

```
library(devtools)
devtools::install_github("NINAnor/oneimpact", ref = "HEAD")
```

## Run with Docker

```bash
docker run --rm -p 8787:8787 -e PASSWORD=rstudio -v $PWD/myproject:/home/rstudio/myproject ghcr.io/ninanor/oneimpact:main
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
  - Get citation information for `oneimpact` in R running `citation(package = 'oneimpact')`
  - Contributions are mostly welcome!
