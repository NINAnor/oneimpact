% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{sample_area_roads.tif}
\alias{sample_area_roads.tif}
\title{Road raster data for the sample area}
\format{
A Geotiff file. Projected CRS: ETRS89 / UTM zone 33N.
\itemize{
\item{1:} {Presence of roads}
\item{NA:} {No presence of roads}
}
}
\source{
\url{https://kartkatalog.geonorge.no/metadata/elveg/ed1e6798-b3cf-48be-aee1-c0d3531da01a}
}
\description{
Raster data indicating pixels with public roads in Southern
Norway. Rasterized from the vector data from the Norwegian road dataset Elveg 1.0
with 100 m resolution and
clipped for the study area presented in the \code{oneimpact} package.
}
\examples{
(f <- system.file("raster/sample_area_roads.tif", package = "oneimpact"))
r <- terra::rast(f)
plot(r)

}
\seealso{
Maps for the sample area: \cr
Limits of sample area: \link{sample_area.gpkg} \cr
Cabins: \link{sample_area_cabins.gpkg}, \link{sample_area_cabins.tif},
\link{sample_area_cabins_count.tif} \cr
Roads: \link{sample_area_roads.gpkg}
}
