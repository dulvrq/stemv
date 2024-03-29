% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stem_volume.R
\name{calcStemVolumeMulti}
\alias{calcStemVolumeMulti}
\title{Calculate stem volume based on junction adjustment, vector version.}
\usage{
calcStemVolumeMulti(
  dt_stem_i,
  D,
  H,
  eq_type = 1,
  adj = "None",
  list_coefs = NULL
)
}
\arguments{
\item{dt_stem_i}{Dataframe from the predefined table.}

\item{D}{Numeric vector. DBH in cm. All elements should be the same DBH class in \code{dt_stem_i}, but it is acceptable
to pass DBH with different classes (processing is a little slow).}

\item{H}{Numeric vector. Tree height in m.}

\item{eq_type}{Numeric. equation type (currently, 1-4 is defined).}

\item{adj}{Adjusting methods for the junction. One of \code{None}, \verb{3w}, \verb{5w}, \verb{3w_o}, \verb{5w_o}.}

\item{list_coefs}{A list of parameters, which is derived from \code{dt_stem_i}. Usually, this is set to NULL.}
}
\value{
Calculated stem volume
}
\description{
Calculate stem volume using an equation corresponding with the equation type.\cr
This function accepts numeric vectors.
}
\details{
This is a wrapper function to apply 3/5- points moving window average for calculating stem volume when
DBH is located at the junction of equations.

See the following reference for calculation details:\cr
Hosoda et al. (2010):
"Differences between the present stem volume tables and the values of the volume equations, and their correction"
Jpn. J. For. Plann. 44:23-39.
}
\references{
\url{https://doi.org/10.20659/jjfp.44.2_23}
}
