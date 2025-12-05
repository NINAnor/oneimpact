# Load a series of files output of fit_net_clogit models and put them on a bag

Load a series of files output of fit_net_clogit models and put them on a
bag

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

  `[logical(1)=TRUE]`  
  If `FALSE` (default), the names of the resamples are taken from the
  `samples` parameter. If `TRUE`, they are taken from the files names,
  instead. In this case, the string defined by the parameter
  `name_from_file_pattern` is used to identify the number of each
  resample.

- name_from_file_pattern:

  `[string="Resample"]`  
  String used to separate the file name and identify the number of the
  resample of each model in the bag, when they are read from files.

  The parameters metric, standardize, and method should be same ones
  used to fit the bag of models.

- f:

  `[formula]`  
  Formula of the models fitted, with all possible candidate terms.

- samples:

  `[list]`  
  List of samples with at least three elements: train, test, and
  validate. Each elements might have several elements, each representing
  the lines of `data` to be sampled for each resample. Typically, this
  is computed by the function
  [`create_resamples()`](https://ninanor.github.io/oneimpact/reference/create_resamples.md).
