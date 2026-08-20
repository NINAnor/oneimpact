# Save kernel/filter matrix to use outside R

This function saves a matrix with weights (filter or kernel matrix) in
an external text file. It can save either the raw matrix or use the
standards for running `r.mfilter` within GRASS GIS (with specific header
and details).

## Usage

``` r
filter_save(
  filt,
  radius,
  type,
  save_format = c("GRASS_rmfilter", "raw")[1],
  save_folder = NULL,
  save_file = NULL,
  divisor = 1,
  normalize = FALSE,
  parallel = TRUE,
  separator = " "
)
```

## Arguments

- filt:

  `[matrix]`  
  Filter or weight matrix, such as one created by
  [`filter_create()`](https://ninanor.github.io/oneimpact/reference/filter_create.md)
  or
  [`terra::focalMat()`](https://rspatial.github.io/terra/reference/focalMat.html).

- radius:

  `[numeric(1)]`  
  Radius of the Zone of Influence (ZOI) of the matrix, in meters.

- type:

  `[character(1)]` Function for the kernel or filter matrix (see `type`
  parameter for
  [`filter_create()`](https://ninanor.github.io/oneimpact/reference/filter_create.md)).

- save_format:

  `[character(1)="GRASS_rmfilter"]{"GRASS_rmfilter", "raw"}`  
  Format in which the function should be saved. Currently, either of the
  two options:

  - GRASS GIS format for the module `r.mfilter`
    (`save_format = "GRASS_rmfilter"`), see details
    [here](https://grass.osgeo.org/grass78/manuals/r.mfilter.html));

  - raw matrix (`save_format = "raw"`), in which only the values of the
    matrix are printed.

- save_folder:

  `[character(1)=NULL]`  
  Path to the folder where the matrix file should be written. If `NULL`,
  the current working directory is used.

- save_file:

  `[character(1)=NULL]`  
  Name of the output file, generally a ".txt" file. If `NULL`, a
  standard filename is created, using the `type` and `radius`. E.g.
  "filter_bartlett2000.txt".

- divisor:

  `[numeric(1)=1]`  
  By default, 1. This is the divisor of the neighborhood matrix when
  used within `r.mfilter`. According the the module documentation, "The
  filter process produces a new category value for each cell in the
  input raster map layer by multiplying the category values of the cells
  in the n x n neighborhood around the center cell by the corresponding
  matrix value and adding them together. If a divisor is specified, the
  sum is divided by this divisor."  
  If the divisor is zero, "then the divisor is computed for each cell as
  the sum of the MATRIX values where the corresponding input cell is
  non-null." In other words, the output map will be rescaled to the
  interval 0,1. If `normalize = TRUE`, the divisor is set to `n*n`.

- normalize:

  `[logical(1)=FALSE]`  
  Whether the matrix should be normalized (sum of all cells is 1 if
  `normalize = TRUE`) or kept as it is (default, `normalize = FALSE`).

- parallel:

  `[logical(1)=TRUE]`  
  Whether the computation should be paralelized or not (details in the
  documentation of the
  [`r.mfilter`](https://grass.osgeo.org/grass78/manuals/r.mfilter.html)
  module).

- separator:

  `[character(1)=" "]`  
  Separator between values of the matrix, within each line. Default is a
  space.

## Value

None. The funcion only saves the input matrix as an external file.

## Details

If used in the `r.mfilter` GRASS GIS module, "The filter process
produces a new category value for each cell in the input raster map
layer by multiplying the category values of the cells in the n x n
neighborhood around the center cell by the corresponding matrix value
and adding them together. If a divisor is specified, the sum is divided
by this divisor." See details
[here](https://grass.osgeo.org/grass78/manuals/r.mfilter.html).

## See also

GRASS module
[r.mfilter](https://grass.osgeo.org/grass78/manuals/r.mfilter.html)

## Examples

``` r
my_filter <- filter_create(r = 100, type = "bartlett", radius = 1000, round = 4)
filter_save(my_filter, radius = 1000, type = "bartlett", save_format = "GRASS_rmfilter")
filter_save(my_filter, radius = 1000, type = "bartlett", save_format = "raw")
```
