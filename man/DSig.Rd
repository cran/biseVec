\name{DSig}
\alias{DSig}

\title{Double Sigmoid Function Creation}
\description{Calculates the green-up day out of true NDVI values extracted from a double sigmoid function fitted to NDVI composites}
\usage{DSig(ndvi, globalthreshold = 0, ravg=FALSE, resvec=c(0.55,0.6,0.65,0.65,0.7))}
\arguments{
  \item{ndvi}{NDVI values computed from the daily surface reflectances, values between 0 and 10000 as a numeric vector.}
  \item{globalthreshold}{The algorithm sets the green-up day at the inflection point of the function if the value is less than 0. 
If the value is 0, the green-up day is calculated through a local threshold method and for a value greater than 0 through a global threshold method.}
  \item{ravg}{If this value is TRUE, a running average will be applied before fitting the double sigmoid function}
  \item{resvec}{...}
}
\details{This routine tries to restore the temporal NDVI profile, i.e. separate true observations from noise. 
An algorithm creates a double sigmoid function fitted to NDVI composites.}
\value{The computed green-up date as Julian day (Day of Year).}
\author{Daniel Doktor, Maximilian Lange}

\examples{
values <- vector(mode="numeric", length=365)
values <- c(seq(from=1000, by=33, length=150), 
	rep(6000, 65), seq(from=6000, by=-33, length=150))
obsday <- 0
obsday <- DSig(values)
obsday
}

