% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{sample_area_cabins.tif}
\alias{sample_area_cabins.tif}
\title{Cabin presence raster data for the sample area}
\format{
A Geotiff file. Projected CRS: ETRS89 / UTM zone 33N.
\itemize{
\item{1:} {Presence of cabins}
\item{NA:} {No presence of cabins}
}
}
\source{
\url{https://register.geonorge.no/det-offentlige-kartgrunnlaget/n50-kartdata/ea192681-d039-42ec-b1bc-f3ce04c189ac}
}
\description{
Raster data indicating pixels with presence of tourist private cabins in Norway.
Cabins corresponds to some specific building types (object_type = "Bygning",
byggtyp_nbr = c("161", "162", "163")) form the public N50 dataset.
The original data consisted of point vector data and were rasterized with 100m
resolution, for the purpose of illustration. The raster was clipped for the
sample area presented in the \code{oneimpact} package.
}
\examples{
(f <- system.file("raster/sample_area_cabins.tif", package = "oneimpact"))
r <- terra::rast(f)
plot(r)

}
\seealso{
Maps for the sample area: \cr
Limits of sample area: \link{sample_area.gpkg} \cr
Cabins: \link{sample_area_cabins.gpkg}, \link{sample_area_cabins_count.tif} \cr
Roads: \link{sample_area_roads.gpkg}, \link{sample_area_roads.tif}
}
