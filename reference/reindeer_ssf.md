# Annotated data of wild reindeer in Norway, prepared for step-selection analysis.

A dataset containing GPS positions for wild reindeer (*Rangifer tarandus
tarandus*) in Setesdal Austhei management area, southern Norway,
organized in a use-availability design for step-selection analysis
(SSA). It includes data from 9 female individuals collected between 2007
and 2010. The coordinates are presented in coordinate reference system
ETRS89/UTM 33N (EPSG:25833). The data set was regularized with 3h data
and 10 random steps were created for each used step. It was then
annotated with rasters on the zone of influence of roads and cabins with
radii 1000 - 3000, both at the start and the end point of a step.

For more information on the original data, see
?[reindeer](https://ninanor.github.io/oneimpact/reference/reindeer.md).

## Usage

``` r
data(reindeer_ssf)
```

## Format

A data frame with 31,735 rows and 28 variables:

- animal_year_id:

  Unique individual-year identifier, used as the sample unit (stratum)
  in the step-selection analysis.

- original_animal_id:

  Animal name assigned at capture.

- sex:

  Sex of the individual animal.

- burst\_:

  Identifier of a continuous tracking burst (uninterrupted sequence of
  GPS fixes for an individual).

- t1\_:

  Timestamp at the start of the step, in UTC.

- t2\_:

  Timestamp at the end of the step, in UTC.

- dt\_:

  Time difference between `t1_` and `t2_`.

- x1\_:

  UTM easting coordinate at the start of the step, in ETRS89/UTM 33N.

- y1\_:

  UTM northing coordinate at the start of the step, in ETRS89/UTM 33N.

- x2\_:

  UTM easting coordinate at the end of the step, in ETRS89/UTM 33N.

- y2\_:

  UTM northing coordinate at the end of the step, in ETRS89/UTM 33N.

- sl\_:

  Step length: Euclidean distance (m) between the start and end of the
  step.

- ta\_:

  Turning angle (radians) between the current and previous step
  direction.

- case\_:

  Case indicator in the use-availability design; 1 = used step, 0 =
  random step.

- step_id\_:

  Step identifier within each burst.

- step_id:

  Global step identifier across all bursts and individuals.

- start_cabins1000:

  Zone of influence of cabins (cumulative) at the start of the step,
  with exponential decay shape and radius 1000 m.

- start_cabins2000:

  Zone of influence of cabins (cumulative) at the start of the step,
  with exponential decay shape and radius 2000 m.

- start_cabins3000:

  Zone of influence of cabins (cumulative) at the start of the step,
  with exponential decay shape and radius 3000 m.

- start_roads1000:

  Zone of influence of roads (cumulative) at the start of the step, with
  exponential decay shape and radius 1000 m.

- start_roads2000:

  Zone of influence of roads (cumulative) at the start of the step, with
  exponential decay shape and radius 2000 m.

- start_roads3000:

  Zone of influence of roads (cumulative) at the start of the step, with
  exponential decay shape and radius 3000 m.

- end_cabins1000:

  Zone of influence of cabins (cumulative) at the end of the step, with
  exponential decay shape and radius 1000 m.

- end_cabins2000:

  Zone of influence of cabins (cumulative) at the end of the step, with
  exponential decay shape and radius 2000 m.

- end_cabins3000:

  Zone of influence of cabins (cumulative) at the end of the step, with
  exponential decay shape and radius 3000 m.

- end_roads1000:

  Zone of influence of roads (cumulative) at the end of the step, with
  exponential decay shape and radius 1000 m.

- end_roads2000:

  Zone of influence of roads (cumulative) at the end of the step, with
  exponential decay shape and radius 2000 m.

- end_roads3000:

  Zone of influence of roads (cumulative) at the end of the step, with
  exponential decay shape and radius 3000 m.

## Source

Panzacchi, M., Van Moorter, B., Strand, O., Loe, L. E., & Reimers, E.
(2015). Searching for the fundamental niche using individual-based
habitat selection modelling across populations. Ecography, 38(7),
659–669. <https://doi.org/10.1111/ecog.01075>

Cagnacci, F., Focardi, S., Ghisla, A., van Moorter, B., Merrill, E.H.,
Gurarie, E., Heurich, M., Mysterud, A., Linnell, J., Panzacchi, M., May,
R., Nygård, T., Rolandsen, C. and Hebblewhite, M. (2016), How many
routes lead to migration? Comparison of methods to assess and
characterize migratory movements. J Anim Ecol, 85: 54-68.
<https://doi.org/10.1111/1365-2656.12449>.

## Examples

``` r
library(tibble)
data("reindeer_ssf")
```
