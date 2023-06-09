% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lazySDM.R
\name{lazySDM}
\alias{lazySDM}
\title{LazySDM}
\usage{
lazySDM(dato, buff = 500, stck1, stck2)
}
\arguments{
\item{dato}{data.frame containing geographic coordinates of a species presence}

\item{buff}{geographic buffer (in km) where to generate pseudo-absence. It is also the area where the model is tested.}

\item{stck1}{Raster stack (terra class) with environmental variables for time 1}

\item{stck2}{Raster stack (terra class) with environmental variables for time 2}
}
\value{
Four objects
summary: Main diagnostic statistics
target.area: raster of the buffer zone
occupancy.time1: raster of estimated habitability for time 1
occupancy.time2: raster of estimated habitability for time 2
}
\description{
Black-box correlational SDM, with minimum specification steps. Calibrates SDM in time1 and projects onto time2.
}
\examples{
lazySDM(whiteshark,buff=100,stck1,stck2)

}
