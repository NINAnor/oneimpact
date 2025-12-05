# Create pretty, standardized sequences from numbers

Create pretty, standardized sequences from numbers

## Usage

``` r
pretty_seq(x, name = "Resample")
```

## Arguments

- x:

  `[vector,numeric]`  
  Vector of integers.

- name:

  `[character(1)="Resample"]`  
  String to be added to the numbers in the sequence.

## Examples

``` r
pretty_seq(1:20)
#>  [1] "Resample01" "Resample02" "Resample03" "Resample04" "Resample05"
#>  [6] "Resample06" "Resample07" "Resample08" "Resample09" "Resample10"
#> [11] "Resample11" "Resample12" "Resample13" "Resample14" "Resample15"
#> [16] "Resample16" "Resample17" "Resample18" "Resample19" "Resample20"
```
