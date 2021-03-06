% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predictionFunctions.R
\name{predCubic}
\alias{predCubic}
\title{predCubic}
\usage{
predCubic(v)
}
\arguments{
\item{v}{is a data frame of values for each paramter in the mcmc temperature model. Required columns are 'Parameter' and 'mean'}
}
\value{
Returns predCubic a data frame with the predictions by day of year and bY which contains the paramter estimates and predicted temperatures from the random effect component of the temperature model
}
\description{
\code{predCubic} Get predictions of the cubic random temperature effect across years
}
\details{
This function extracts the among-year temperature random effects parameters from the model output, does predictions of the effects by day of year and prevodes predictions at 4 day of year values (z-score values)
}
\examples{

\dontrun{
preds <- predCubic( values )
}
}

