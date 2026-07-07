# Load a series of files output of fit_net\_(c)logit models and put them on a bag

This is a helper function aimed at loading a bag of models created
through a call to the
[`fit_net_clogit()`](https://ninanor.github.io/oneimpact/reference/fit_net_functions.md)
or
[`bag_fit_net_clogit()`](https://ninanor.github.io/oneimpact/reference/bag_fit_net_functions.md)
functions, when the output is saved in external files. Files saved in
`.rds` format are here loaded together so they can be provided to the
[`bag_models()`](https://ninanor.github.io/oneimpact/reference/bag_models.md)
function to build a `bag` object with the fitted models.

## Usage

``` r
bag_load_models(
  data,
  load_models_path = ".",
  load_models_pattern = NULL,
  names_from_file = FALSE,
  name_from_file_pattern = "Resample",
  verbose = FALSE
)
```

## Arguments

- data:

  `[data.frame,tibble]`  
  Complete data set analyzed.

- load_models_path:

  `[string="."]`  
  Path to the folder where the files are saved.

- load_models_pattern:

  `[string="."]`  
  Pattern common to the file names, to be used in
  [`grep()`](https://rdrr.io/r/base/grep.html) to find the files within
  the folder `load_models_path`. It should be the pattern present in the
  argument `out_dir_file` in
  [`fit_net_clogit()`](https://ninanor.github.io/oneimpact/reference/fit_net_functions.md)
  or
  [`bag_fit_net_clogit()`](https://ninanor.github.io/oneimpact/reference/bag_fit_net_functions.md)
  functions.

- names_from_file:

  `[logical(1)=FALSE]`  
  If `FALSE` (default), the names of the resamples are taken from the
  loaded model objects. If `TRUE`, they are taken from the files names,
  instead. In this case, the string defined by the parameter
  `name_from_file_pattern` is used to identify the number of each
  resample.

- name_from_file_pattern:

  `[string="Resample"]`  
  String used to separate the file name and identify the number of the
  resample of each model in the bag, when they are read from files.

- verbose:

  `[logical(1)=FALSE]`  
  Should messages of the computation steps be printed in the prompt
  along the computation?

## Value

A `list` with loaded model objects and metadata for bag construction, in
the same format as if all the models were fitted by the
[`fit_net_clogit()`](https://ninanor.github.io/oneimpact/reference/fit_net_functions.md)
or
[`bag_fit_net_clogit()`](https://ninanor.github.io/oneimpact/reference/bag_fit_net_functions.md)
functions within a R session. Key fields include:

- `n`: expected number of resamples/models.

- `formula`: model formula used in the loaded fits.

- `method`: fitting method used (e.g., `"Lasso"`, `"AdaptiveLasso"`).

- `metric`: validation metric used in the loaded fits.

- `samples`: resampling structure from the fitted objects.

- `standardize`: standardization mode used during fitting.

- `models`: named list of loaded model objects (from `.rds` files).

This object is intended as input to
[`bag_models()`](https://ninanor.github.io/oneimpact/reference/bag_models.md).
