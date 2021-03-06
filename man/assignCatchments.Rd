<<<<<<< HEAD
% Generated by roxygen2 (4.0.2): do not edit by hand
=======
% Generated by roxygen2 (4.0.1): do not edit by hand
>>>>>>> 34a3efc373678e8b541effdb2d515ef4aef95807
\name{assignCatchments}
\alias{assignCatchments}
\title{assignCatchments}
\usage{
assignCatchments(sites, catchmentShapefile, catchmentID, projectionString)
}
\arguments{
\item{sites}{dataframe containing 3 columns in specified order: 1)SiteID 2)Longitude 3)Latitude}

\item{catchmentShapefile}{SpatialPolygonsDataframe of the catchment shapefile over which the points will be matched}

\item{catchmentID}{a character vector of column name describing the catchment identifier in the shapefile (e.g. "FEATUREID" for NHDplusV2)}

\item{projectionString}{a CRS string of the spatial data projection of the shapefile and site coordinates}
}
\value{
Returns the input dataframe with the catchment ID appended (column 4)
}
\description{
\code{assignCatchments} Assigns a catchment to a site with a set of lat/lon points.
}
\details{
This function uses a spatial overlay to assign catchment IDs to sites with associated lat/lon points.
If a point does not match any catchment, NA is returned
}

