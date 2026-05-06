# Package index

## Zone of influence (ZoI) functions

Functions to represent the zone of influence

- [`dist_decay()`](https://ninanor.github.io/oneimpact/reference/zoi_functions.md)
  [`threshold_decay()`](https://ninanor.github.io/oneimpact/reference/zoi_functions.md)
  [`step_decay()`](https://ninanor.github.io/oneimpact/reference/zoi_functions.md)
  [`bartlett_decay()`](https://ninanor.github.io/oneimpact/reference/zoi_functions.md)
  [`tent_decay()`](https://ninanor.github.io/oneimpact/reference/zoi_functions.md)
  [`linear_decay()`](https://ninanor.github.io/oneimpact/reference/zoi_functions.md)
  [`gaussian_decay()`](https://ninanor.github.io/oneimpact/reference/zoi_functions.md)
  [`half_norm_decay()`](https://ninanor.github.io/oneimpact/reference/zoi_functions.md)
  [`exp_decay()`](https://ninanor.github.io/oneimpact/reference/zoi_functions.md)
  : Zone of Influence (ZOI) functions
- [`plot_zoi1d()`](https://ninanor.github.io/oneimpact/reference/plot_zoi1d.md)
  : Plot the functions for the nearest and cumulative zone of influence
  in 1 dimension

## Compute zones of influence (ZoI)

Functions to compute the zone of influence for raster maps

- [`calc_zoi()`](https://ninanor.github.io/oneimpact/reference/calc_zoi.md)
  : Calculates the zone of influence from the nearest feature and the
  cumulative zone of influence of multiple features
- [`calc_zoi_cumulative()`](https://ninanor.github.io/oneimpact/reference/calc_zoi_cumulative.md)
  : Calculate the cumulative zone of influence of multiple features
- [`calc_zoi_nearest()`](https://ninanor.github.io/oneimpact/reference/calc_zoi_nearest.md)
  : Calculate the zone of influence from the nearest feature
- [`calc_zoi_sql()`](https://ninanor.github.io/oneimpact/reference/calc_zoi_sql.md)
  : Compute Zone of Influence for points and annotate them to points
  using SQL

## Spatial filters

Creating spatial filters to compute a meaningful zone of influence

- [`filter_create()`](https://ninanor.github.io/oneimpact/reference/filter_create.md)
  : Create filters or kernel matrices for raster neighborhood analyses
- [`filter_na_strata()`](https://ninanor.github.io/oneimpact/reference/filter_na_strata.md)
  : Remove missing values and ensure case-control per stratum
- [`filter_save()`](https://ninanor.github.io/oneimpact/reference/filter_save.md)
  : Save kernel/filter matrix to use outside R

## Estimating ZOI - set up analysis

Functions to set up RSF and SSF analyses using ZOI variables

- [`add_zoi_formula()`](https://ninanor.github.io/oneimpact/reference/add_zoi_formula.md)
  : Adds ZOI radii to formula
- [`add_zoi_layers()`](https://ninanor.github.io/oneimpact/reference/add_zoi_layers.md)
  : Create ZOI layer names as strings for data annotation
- [`explore_blocks()`](https://ninanor.github.io/oneimpact/reference/explore_blocks.md)
  : Explore hierarchical blocks after sampling or spatial stratification
- [`explore_blocks_pre()`](https://ninanor.github.io/oneimpact/reference/explore_blocks_pre.md)
  : Explore potential hierarchical blocks before sampling or spatial
  stratification
- [`spat_strat()`](https://ninanor.github.io/oneimpact/reference/spat_strat.md)
  : Preparing data for spatially stratified cross‐validation schemes
- [`create_resamples()`](https://ninanor.github.io/oneimpact/reference/create_resamples.md)
  : Create samples for fitting, calibrating, and validating models

## Estimating ZOI - fit models

Functions to fit RSF and SSF and estimate ZOI using penalized regression

- [`filter_na_strata()`](https://ninanor.github.io/oneimpact/reference/filter_na_strata.md)
  : Remove missing values and ensure case-control per stratum
- [`extract_response_strata()`](https://ninanor.github.io/oneimpact/reference/extract_response_strata.md)
  : Separates elements in a statistical formula
- [`net_clogit()`](https://ninanor.github.io/oneimpact/reference/net_clogit.md)
  [`net_ssf()`](https://ninanor.github.io/oneimpact/reference/net_clogit.md)
  [`net_issf()`](https://ninanor.github.io/oneimpact/reference/net_clogit.md)
  : Fits a conditional logistic regression/SSF/iSSF using glmnet
- [`net_logit()`](https://ninanor.github.io/oneimpact/reference/net_logit.md)
  [`net_rsf()`](https://ninanor.github.io/oneimpact/reference/net_logit.md)
  : Fits a logistic regression/RSF using glmnet
- [`fit_net_clogit()`](https://ninanor.github.io/oneimpact/reference/fit_net_functions.md)
  [`fit_net_ssf()`](https://ninanor.github.io/oneimpact/reference/fit_net_functions.md)
  [`fit_net_issf()`](https://ninanor.github.io/oneimpact/reference/fit_net_functions.md)
  [`fit_net_logit()`](https://ninanor.github.io/oneimpact/reference/fit_net_functions.md)
  [`fit_net_rsf()`](https://ninanor.github.io/oneimpact/reference/fit_net_functions.md)
  [`grouped_func()`](https://ninanor.github.io/oneimpact/reference/fit_net_functions.md)
  : Fits a conditional logistic regression/SSF/iSSF with penalized
  regression using glmnet in a train-validate-test setup
- [`bag_fit_net_clogit()`](https://ninanor.github.io/oneimpact/reference/bag_fit_net_functions.md)
  [`bag_fit_net_logit()`](https://ninanor.github.io/oneimpact/reference/bag_fit_net_functions.md)
  : Fit a bag of conditional logistic regression/SSF/iSSF models with
  penalized regression in a train-validate-test setup
- [`bag_load_models()`](https://ninanor.github.io/oneimpact/reference/bag_load_models.md)
  : Load a series of files output of fit_net_clogit models and put them
  on a bag
- [`bag_models()`](https://ninanor.github.io/oneimpact/reference/bag_models.md)
  : Summary of a bag of models
- [`truncate_bag()`](https://ninanor.github.io/oneimpact/reference/truncate_bag.md)
  : Truncate bag to avoid weirdness in the model
- [`conditionalBoyce()`](https://ninanor.github.io/oneimpact/reference/concordance_indices.md)
  [`conditionalSomersD()`](https://ninanor.github.io/oneimpact/reference/concordance_indices.md)
  [`conditionalAUC()`](https://ninanor.github.io/oneimpact/reference/concordance_indices.md)
  [`AUC()`](https://ninanor.github.io/oneimpact/reference/concordance_indices.md)
  [`coxnet.deviance()`](https://ninanor.github.io/oneimpact/reference/concordance_indices.md)
  [`Cindex()`](https://ninanor.github.io/oneimpact/reference/concordance_indices.md)
  : Computes the conditional Boyce index for model evaluation
- [`kernel_prediction()`](https://ninanor.github.io/oneimpact/reference/kernel_prediction.md)
  : Prediction based only on step length terms

## Estimating ZOI - interpret and visualize models

Functions to help interpreting parameters and visualizing cumulative
impacts

- [`predict()`](https://ninanor.github.io/oneimpact/reference/predict.md)
  : Prediction of a bag of models to new data
- [`zoi_from_curve()`](https://ninanor.github.io/oneimpact/reference/zoi_from_curve.md)
  : Get estimates of zone of influence (ZOI) metrics from response
  curves
- [`plot_importance()`](https://ninanor.github.io/oneimpact/reference/plot_importance.md)
  : Plot variable importance
- [`variable_importance()`](https://ninanor.github.io/oneimpact/reference/variable_importance.md)
  : Computes variable importance
- [`plot_coef()`](https://ninanor.github.io/oneimpact/reference/plot_coef.md)
  : Plot the coefficients of bags of models
- [`plot_response()`](https://ninanor.github.io/oneimpact/reference/plot_response.md)
  : Plot responses from a bag of models
- [`bag_predict_spat()`](https://ninanor.github.io/oneimpact/reference/bag_predict_spat.md)
  [`bag_predict_spat_vars()`](https://ninanor.github.io/oneimpact/reference/bag_predict_spat.md)
  : Predict bag of models in space
- [`implausibility()`](https://ninanor.github.io/oneimpact/reference/implausibility.md)
  : Computes ecological implausibility for a fitted model or its
  estimated coefficients
- [`rescale_coefficients()`](https://ninanor.github.io/oneimpact/reference/rescale_coefficients.md)
  : Rescale standardized coefficients back to their original range after
  model fitting

## Raster processing ancillary functions

GRASS GIS ancillary functions to process rasters and prepare inputs for
computing zones of influence

- [`grass_binarize()`](https://ninanor.github.io/oneimpact/reference/grass_binarize.md)
  : Binarize continuous raster maps
- [`grass_find_layer()`](https://ninanor.github.io/oneimpact/reference/grass_find_layer.md)
  : Find layers within GRASS GIS with multiple patterns
- [`grass_v2rast_count()`](https://ninanor.github.io/oneimpact/reference/grass_v2rast_count.md)
  : Rasterizes a vector counting the number of features in each pixel
- [`raster_rescale()`](https://ninanor.github.io/oneimpact/reference/raster_rescale.md)
  : Rescale raster values

## Support for landscape simulation

- [`set_points()`](https://ninanor.github.io/oneimpact/reference/set_points.md)
  : Simulate points in a landscape
- [`set_points_from_raster()`](https://ninanor.github.io/oneimpact/reference/set_points_from_raster.md)
  : Simulate points using input raster as weights
- [`set_points_sample()`](https://ninanor.github.io/oneimpact/reference/set_points_sample.md)
  : Simulate regular or random points in 2D
- [`isolation()`](https://ninanor.github.io/oneimpact/reference/isolation.md)
  [`mean_isolation()`](https://ninanor.github.io/oneimpact/reference/isolation.md)
  : Isolation and mean isolation of points in space

## Datasets

Datasets for testing the ZoI approach

- [`reindeer`](https://ninanor.github.io/oneimpact/reference/reindeer.md)
  : GPS positions for wild reindeer in Norway.
- [`reindeer_area.gpkg`](https://ninanor.github.io/oneimpact/reference/reindeer_area.gpkg.md)
  : Reindeer area: a polygon vector data for the Setesdal Austhei
  reindeer herding area
- [`reindeer_cabins.gpkg`](https://ninanor.github.io/oneimpact/reference/reindeer_cabins.gpkg.md)
  : Cabins vector data for the reindeer area
- [`reindeer_roads_private.gpkg`](https://ninanor.github.io/oneimpact/reference/reindeer_roads_private.gpkg.md)
  : Private roads vector data for the reindeer area
- [`reindeer_roads_public.gpkg`](https://ninanor.github.io/oneimpact/reference/reindeer_roads_public.gpkg.md)
  : Public roads vector data for the reindeer area
- [`reindeer_rsf`](https://ninanor.github.io/oneimpact/reference/reindeer_rsf.md)
  : Annotated data of wild reindeer in Norway, prepared for point
  resource-selection analysis.
- [`reindeer_ssf`](https://ninanor.github.io/oneimpact/reference/reindeer_ssf.md)
  : Annotated data of wild reindeer in Norway, prepared for
  step-selection analysis.
- [`rast_predictors_hardanger_500m.tif`](https://ninanor.github.io/oneimpact/reference/rast_predictors_hardanger_500m.tif.md)
  : Predictor environmental variables for the Hardangervidda wild
  reindeer area in Norway
- [`sample_area_cabins.tif`](https://ninanor.github.io/oneimpact/reference/sample_area_cabins.tif.md)
  : Cabin presence raster data for the sample area
- [`sample_area_cabins_count.tif`](https://ninanor.github.io/oneimpact/reference/sample_area_cabins_count.tif.md)
  : Cabin count raster data for the sample area
- [`sample_area_roads.tif`](https://ninanor.github.io/oneimpact/reference/sample_area_roads.tif.md)
  : Road raster data for the sample area
- [`sample_area.gpkg`](https://ninanor.github.io/oneimpact/reference/sample_area.gpkg.md)
  : Sample area: a polygon vector data
- [`sample_area_cabins.gpkg`](https://ninanor.github.io/oneimpact/reference/sample_area_cabins.gpkg.md)
  : Cabins vector data for the sample area
- [`sample_area_roads.gpkg`](https://ninanor.github.io/oneimpact/reference/sample_area_roads.gpkg.md)
  : Road vector data for the sample area
