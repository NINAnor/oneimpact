# read data
data("reindeer_ssf")

# explore blocks - animal ID as block H0
explore_blocks_pre(reindeer_ssf, "original_animal_id", col_case = "case_")

# explore blocks - year as block H0
library(lubridate)
reindeer_ssf |>
  dplyr::mutate(year = lubridate::year(t1_)) |>
  explore_blocks_pre("year", col_case = "case_")

# year as block H0 + animal ID
reindeer_ssf |>
  dplyr::mutate(year = lubridate::year(t1_)) |>
  explore_blocks_pre("year", animal_id = "original_animal_id", col_case = "case_")
