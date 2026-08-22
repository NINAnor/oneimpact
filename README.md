# oneimpact <img src="man/figures/oneimpact_hex_logo.png" align="right" alt="" width="150" />

<!-- badges: start -->
  [![DOI](https://zenodo.org/badge/453101311.svg)](https://zenodo.org/badge/latestdoi/453101311)
<!--  [![R-CMD-check](https://github.com/NINAnor/oneimpact/workflows/R-CMD-check/badge.svg)](https://github.com/NINAnor/oneimpact/actions) -->
<!-- badges: end -->

`oneimpact` provides tools for the assessment of cumulative impacts of multiple infrastructure and land use modifications in ecological studies. This includes tools to calculate the zone of influence (ZOI, also called scale of effect in landscape ecology) of anthropogenic variables as well as tools for model fitting, estimation of the effect size and ZOI, making predictions and interpreting fitted relationships, computing model performance and plausibility, and ancillary functions. The functions dealing with spatial data processing can be run in both R and GRASS GIS, using R as an interface.

## Context: Why zone of influence and cumulative impacts?

Traditional ecological impact assessments typically focus on the effects of **distance to the nearest feature only**, which might underestimate impacts where multiple anthropogenic stressors are clustered in space. 
The **cumulative ZOI approach** implemented in `oneimpact` recognizes that:

1. **The impact is a product of both effect size and zone of influence**: impacts extend beyond the features themselves — a road, cabin, or wind farm may affect species and processes across the landscapes, not just at its location.
2. **The impact of multiple features might accumulate**: impacts of features clustered in space (e.g., "cabin villages" vs. single, isolated cabins) can exceed impacts of isolated features in magnitude and spatial extent.
3. **the ZOI is context-dependent**: The spatial extent of an impact varies by species, process, and disturbance type.

For detailed background and comprehensive workflow examples, see the package vignettes.

## Scientific foundation

This package implements frameworks from two peer-reviewed studies:

- **Niebuhr et al. (2023)**: [Estimating the cumulative impact and zone of influence of anthropogenic features on biodiversity](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13994). *Methods in Ecology and Evolution*, 14(6), 1470–1487.
  - Introduces the cumulative ZOI approach with mathematical framework and empirical case study

- **Niebuhr et al. (2026)**: Ecology-informed machine learning to estimate the zone of influence of multiple disturbances on ecological niche models. *Working manuscript*.
  - Extends framework with penalized regression, ecological constraints, and bootstrap aggregation for fitting complex datasets

## Typical use cases

- Habitat suitability modeling: how human infrastructure affects species persistence
- Landscape permeability: assessing movement probability fragmented by anthropogenic disturbance
- Species distribution modeling: Estimating responses to cumulative anthropogenic pressures
- Environmental impact assessment: Predicting magnitude (effect size) and extent (ZOI) of proposed development impacts
- Conservation planning: Identifying priority areas based on cumulative impacts

## Quick Start

```R
# library(devtools)
# devtools::install_github("NINAnor/oneimpact", ref = "HEAD")
library(oneimpact)

# 1. Compute ZOI layers (nearest and cumulative)
zoi <- calc_zoi(
  x = infrastructure_raster,
  radius = c(500, 1000, 2000),
  type = c("linear", "gaussian")
)

# 2. Prepare data with ZOI terms
formula <- add_zoi_formula(
  f = ~slope + vegetation + infrastructureXXX,        # base formula
  zoi_radius = c(500, 1000, 2000),
  type = c("linear", "gaussian"),
  pattern = "XXX",                # placeholder in raster names
  predictor_table = TRUE          # returns both formula and variable table
)$formula

# 3. Set up resampling scheme, random cross-validation
samples <- create_resamples(
  y = your_data$response,
  p = c(0.2, 0.2, 0.2),          # train, validate, test fractions
  times = 100                    # number of resamples
)

# 4. Fit ensemble model
fitted_models <- bag_fit_net_logit(
  f = formula,
  data = your_data,
  samples = samples,
  metric = "AUC",
  method = "Lasso",
  standardize = "internal"
)

# 5. Create bag object
bag <- bag_models(fitted_models, your_data, score_threshold = 0.7)

# 6. Evaluate and visualize
plot_response(bag, predictor = "zoi_500_linear")
variable_importance(bag)
implausibility(bag, response_variable = your_data$response)
impact_map <- bag_predict_spat(bag, newdata = raster_stack)
```

For complete workflow examples, see the vignettes.

## Compute spatial layers representing zones of influence

The first set of functions in `oneimpact` are aimed at computing the (potential) ZOI of infrastructure or other
spatial covariates. This means we use spatial information on where they are located to compute
the density of features in space (i.e. the cumulative ZOI) and/or the (decay) distance to the nearest
feature (i.e. the ZOI of the nearest), given an expected ZOI radius (i.e. the distance up to which
a given feature is expected to affect a certain species or process). These functions do not estimate
the ZOI, though (which is context and process dependent); for that see more functions further down.

Here are the main functions in `oneimpact` to compute spatial layers representing zones of influence.

### Zone of influence (ZOI) decay functions

- zoi_functions: a set of decay zone of influence functions to characterize different shapes of the ZOI around infrastructure, 
parameterized based on the zone of influence radius. The functions implemented so far are: threshold (`threshold_decay()` or `step_decay()`),
linear decay (`linear_decay()` or `bartlett_decay()` or `tent_decay()`), exponential decay (`exp_decay()`), or Gaussian decay 
(`gaussian_decay()` or `half_norm_decay()`).
- `plot_zoi1d()`: plot ZOI in 1 dimensional space for multiple points infrastructure, using both the ZOI of the nearest
feature and the cumulative ZOI metric.

### Compute zones of influence (ZOI)

- `calc_zoi_nearest()`: Calculate the zone of influence from the nearest infrastructure, according to multiple possible 
decay functions and zones of influence radii.
- `calc_zoi_cumulative()`: Calculate the cumulative zone of influence of multiple features, according to multiple possible 
decay functions and zones of influence radii.
- `calc_zoi()`: Calculate both the the ZOI of the nearest infrastructure and the cumulative ZOI, at multiple
scales or zones of influence radii.

### Spatial filters

- `create_filter()`: Create filters or weight matrices for neighborhood analysis, according to different decay functions
and parameterized using the zone of influence radius.
- `save_filter()`: Saves filters/weight matrices outside R for use within GRASS GIS modules.

## Estimate the cumulative impact and the ZOI of features on a certain species or process

The `oneimpact` package also allows us to, given a set of potential candidate ZOIs (with possibly
different types, shapes, and radii; computed with the functions above), estimate the actual 
effect and ZOI of the variables on a certain species or process. This is done combining three elements:

- Bootstrap aggregation (bagging), a multi-model bootstrap procedure that allows us to estimate the 
uncertainty in the effect sizes and ZOI radii;
- Penalized regression, an approach that allows us to penalize estimated coefficients and possibly 
remove the least likely covariates from a model, i.e., it allows us to perform model fitting together
with variable section;
- Nested cross-validation, which allows is to consider hierarchical, spatial, or temporally cross-validation
schemes in model and variable/feature selection.

### Estimating ZOI - set up analysis

Functions to set up RSF and SSF analyses using ZOI variables:

- `add_zoi_formula()`: Adds ZOI radii to formula
- `spat_strat()`: Prepares data for spatially stratified cross‐validation schemes
- `explore_blocks_pre()` and `explore_blocks()`: Explore hierarchical blocks before or 
after sampling or spatial stratification, respectively
- `create_resamples()`: Create samples for fitting, calibrating, and validating models in
a bootstrap/baggin procedure.

### Estimating ZOI - fit models

Functions to fit RSF and SSF and estimate ZOI using penalized regression

- `bag_fit_net_clogit()`: Fits a a bag of conditional logistic regressions/SSF/iSSF using glmnet. 
This function is a wrapper around `fit_net_clogit()` which is the one properly setting up the 
model fitting, tunning, and validation. It allows the use of different penalization algorithms,
including Lasso, Ridge, Adaptive Lasso, and different adaptations from Adaptive Lasso. This
function calls the function `net_logit()` which is the one properly calling `glmnet` and fitting
the model.
- `bag_fit_net_logit()` (and `fit_net_logit()`, `net_logit()`): equivalent to the one above, but
performing common logistic regression, with no strata.

- `bag_load_models()`: Load a vector of files with the output of `fit_net_clogit()` or 
`fit_net_logit()` and put them on a bag.
- `bag_models()`: Bag a list of loaded/fitted models fitted through `fit_net_clogit()` or 
`fit_net_logit()`. This created an object of class `bag` with all information for understanding
and making prediction from the bag of models.
- `AUC()`, `conditionalAUC()`, `coxnet.deviance()`, `Cindex()`, `conditionalSomersD()`: functions
used for model tunning (selecting penalties) and validation.

### Estimating ZOI - interpret and visualize models

Functions to help interpreting parameters and visualizing cumulative impacts from bags
of fitted models:

- `predict()`: Prediction of a bag of models to new data.
- `variable_importance()`, `plot_importance()`: Computes and plots variable importance from a bag of models.
- `plot_coef()`: Plots the coefficients of bags of models.
- `plot_response()`: Plots (partial) response curves from a bag of models.
- `bag_predict_spat()`: Predict bag of models in space.
- `bag_predict_spat_vars()`: Predict reponses of each individual covariate in space according to 
a bag of models.

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

If you use Compose:

```bash
docker compose run rstudio
```

You can customize `docker-compose.yml` based on your needs.

## See also

For model fitting and estimation of ZOI, see the pacakage [`glmnet`](https://glmnet.stanford.edu/index.html),
which is the backbone of the modeling approach used in `oneimpact`. For other similar approaches,
check the [`maxnet()`](https://cran.r-project.org/web/packages/maxnet/index.html) for MaxEnt 
species distribution models using `glmnet`.

The `oneimpact` functions to compute the ZOI layers are greatly based on neighborhood analyses 
made through the [`terra` package](https://rspatial.org/terra/pkg/index.html) in R and on three GRASS GIS modules:
[`r.mfilter`](https://grass.osgeo.org/grass78/manuals/r.mfilter.html), 
[`r.resamp.filter`](https://grass.osgeo.org/grass78/manuals/r.resamp.filter.html), and 
[`r.neighbors`](https://grass.osgeo.org/grass78/manuals/r.neighbors.html). The connection
between R and GRASS GIS is made through the [`rgrass`](https://github.com/rsbivand/rgrass) R package.

## Meta

  - Please [report any issues or bugs](https://github.com/NINAnor/oneimpact/issues/new/).
  - License: GPL3
  - Get citation information for `oneimpact` in R running `citation(package = 'oneimpact')`, or check the reference [here](https://ninanor.github.io/oneimpact/authors.html#citation).
  - Contributions are mostly welcome!
