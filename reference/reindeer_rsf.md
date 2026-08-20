# Annotated data of wild reindeer in Norway, prepared for point resource-selection analysis.

A data set where GPS positions of wild reindeer (*Rangifer tarandus
tarandus*) Hardangervidda management area, southern Norway, were
organized in a use-availability design for point resource-selection
analysis (RSA) and annotated with environmental data. It includes data
from 115 female individuals collected between 2001 and 2019. The data
set was regularized with a 3h fix rate and 9 random locations were
created for each used location. It was then annotated with rasters on
the zone of influence of private cabins, public resorts, and tourist
trails, with exponential decay shape and radii varying from 100m to 10
km, as well as land cover and bio-geo-climatic variables. Geographical
coordinates of the used and random positions were omitted after data
annotation.

This is part of the data set used for analysis in Niebuhr et al. (2023).
ZOI variables with shapes other than the exponential decay were omitted.

## Usage

``` r
data(reindeer_rsf)
```

## Format

A data frame with 74,337 rows and 88 variables:

- use:

  Case in the use-availability setup; 1 represents a used location, 0
  represents a random locations.

- norway_pca_klima_axis1-4:

  Components 1 to 4 from a principal component analysis representing
  bio-geo-climatic variation in Norway, from Bakkestuen et al. 2008.
  PCAs 1 to 4 represent, respectively, continentality, altitude, terrain
  ruggedness, and solar radiation. More information in Niebuhr et al.
  2023.

- norway_pca_klima_axis1-2_sq:

  Squared value for components 1 and 2 to from a principal component
  analysis representing bio-geo-climatic variation in Norway, from
  Bakkestuen et al. 2008.

- NORUTreclass:

  Land use and land cover classes from NORUT, reclassified as in Niebuhr
  et al. 2023.

- cabins_private_nearest_exp_decayXXX:

  Zone of influence of the nearest private cabin at each location, with
  exponential decay shape, and radii defined by XXX (from 100 to 10000
  m).

- cabins_private_cumulative_exp_decayXXX:

  Cumulative zone of influence of private cabins at each location, with
  exponential decay shape, and radii defined by XXX (from 100 to 10000
  m).

- cabins_public_nearest_exp_decayXXX:

  Zone of influence of the nearest public resort at each location, with
  exponential decay shape, and radii defined by XXX (from 100 to 10000
  m).

- cabins_public_cumulative_exp_decayXXX:

  Cumulative zone of influence of public resorts at each location, with
  exponential decay shape, and radii defined by XXX (from 100 to 10000
  m).

- trails_nearest_exp_decayXXX:

  Zone of influence of the nearest trail, weighted by the number of
  tourists in the trail, with exponential decay shape, and radii defined
  by XXX (from 100 to 10000 m).

- trails_cumulative_exp_decayXXX:

  Cumulative zone of influence of trails, weighted by the number of
  tourists in each trail, with exponential decay shape, and radii
  defined by XXX (from 100 to 10000 m).

## Source

Niebuhr, B. B., Van Moorter, B., Stien, A., Tveraa, T., Strand, O.,
Langeland, K., Sandström, P., Alam, M., Skarin, A., & Panzacchi, M.
(2023). Estimating the cumulative impact and zone of influence of
anthropogenic features on biodiversity. Methods in Ecology and
Evolution. https://doi.org/10.1111/2041-210X.14133

Bakkestuen, V., Erikstad, L., & Halvorsen, R. (2008). Step-less models
for regional environmental variation in Norway. Journal of Biogeography,
35(10), 1906–1922. https://doi.org/10.1111/j.1365-2699.2008.01941.x

## Examples

``` r
library(tibble)
data("reindeer_rsf")
```
