% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets-docs.R
\docType{data}
\encoding{UTF-8}
\name{fdir_sub_areas}
\alias{fdir_sub_areas}
\title{Norwegian sub-areas (lokasjon) for commercial fishing}
\format{
\code{\link[sf:st_sf]{sf object}} containing major fishing zones defined by the Norwegian Directorate of Fisheries.
}
\source{
\href{https://open-data-fiskeridirektoratet-fiskeridir.hub.arcgis.com/}{Norwegian Directorate of Fisheries}
}
\usage{
fdir_sub_areas
}
\description{
Norwegian sub-areas (lokasjon) for commercial fishing
}
\examples{
if(requireNamespace("ggspatial")) {
\donttest{
basemap(fdir_sub_areas) + 
ggspatial::annotation_spatial(fdir_sub_areas, fill = NA)
}
}
}
\seealso{
Other datasets: 
\code{\link{fdir_main_areas}},
\code{\link{ices_areas}}
}
\concept{datasets}
\keyword{datasets}
\keyword{fishery}
\keyword{shapefiles}
