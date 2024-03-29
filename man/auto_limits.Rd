% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auto_limits.R
\name{auto_limits}
\alias{auto_limits}
\title{Automatic limits for basemap}
\usage{
auto_limits(
  data,
  lon = NULL,
  lat = NULL,
  proj.in = 4326,
  proj.out = NULL,
  expand.factor = NULL,
  verbose = FALSE
)
}
\arguments{
\item{data}{Data frame or a spatial object containing data for which the limits should be calculated.}

\item{lon, lat}{Names of longitude and latitude columns in \code{data} as character or integer index. If \code{NULL}, the column names are \link[=guess_coordinate_columns]{guessed}.}

\item{proj.in}{Original \code{\link[sf:st_crs]{CRS}} projection. Must be defined as character argument.}

\item{proj.out}{Resulting map projection. See \code{\link{transform_coord}}.}

\item{expand.factor}{Expansion factor for map limits. Set to \code{NULL} to ignore.}

\item{verbose}{Logical indicating whether information about the projection and guessed column names should be returned as message. Set to \code{FALSE} to make the function silent.}
}
\value{
A list of limits and projections in \code{proj.in} and \code{proj.out} formats.
}
\description{
Find limits for a \code{\link{basemap}} from a data frame.
}
\details{
This is an internal function, which is automatically run by the \code{\link{basemap}} function.
}
\examples{
auto_limits(data = expand.grid(lon = c(-120, 180, 120),
   lat = c(60, 60, 80)))
}
\seealso{
Other customize shapefiles: 
\code{\link{reorder_layers}()},
\code{\link{theme_map}()}
}
\author{
Mikko Vihtakari
}
\concept{customize shapefiles}
\keyword{internal}
