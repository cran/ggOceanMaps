% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets-docs.R
\docType{data}
\encoding{UTF-8}
\name{fdir_main_areas}
\alias{fdir_main_areas}
\title{Major fisheries areas (hovedomraade) of Norway}
\format{
\code{\link[sf:st_sf]{sf object}} containing major fishing zones defined by the Norwegian Directorate of Fisheries. Contains also Northwest Atlantic Fisheries Organization's divisions where Norwegian vessels tend to fish.
}
\source{
\href{https://open-data-fiskeridirektoratet-fiskeridir.hub.arcgis.com/}{Norwegian Directorate of Fisheries} and \href{https://www.nafo.int/About-us/Maps}{Northwest Atlantic Fisheries Organization}
}
\usage{
fdir_main_areas
}
\description{
Major fisheries areas (hovedomraade) of Norway
}
\examples{
if(requireNamespace("ggspatial")) {
\donttest{ 
basemap(fdir_main_areas) + 
ggspatial::annotation_spatial(fdir_main_areas, fill = NA)
}
}
}
\seealso{
Other datasets: 
\code{\link{fdir_sub_areas}},
\code{\link{ices_areas}}
}
\concept{datasets}
\keyword{datasets}
\keyword{fishery}
\keyword{shapefiles}
