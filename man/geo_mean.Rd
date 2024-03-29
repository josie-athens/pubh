% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/descriptive_functions.R
\name{geo_mean}
\alias{geo_mean}
\title{Geometric mean.}
\usage{
geo_mean(x)
}
\arguments{
\item{x}{A numeric variable with no negative values.}
}
\value{
A scalar, the calculated geometric mean.
}
\description{
Geometric mean.
}
\examples{
data(IgM, package = "ISwR")
Ab <- data.frame(IgM)
estat(~IgM, data = Ab)
geo_mean(IgM)
}
