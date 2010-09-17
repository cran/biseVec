\name{gauss}
\alias{gauss}

\title{Gaussian Function Creation}
\description{Calculates the green-up day out of true NDVI values extracted from a gaussian function fitted to NDVI composites}
\usage{gauss(ndvi, globalthreshold = 0, asym=FALSE, resvec=c(0.55,0.6,0.65,0.65,0.7))}
\arguments{
  \item{ndvi}{NDVI values computed from the daily surface reflectances, values between 0 and 10000 as a numeric vector.}
  \item{globalthreshold}{The algorithm returns the green-up date calculated through a local threshold method if the value of mode is 0 or through a global threshold method (with a global threshold equal to the value of globalthreshold) if the value is greater than 0.}
  \item{asym}{An asymmetric Gaussian function will be created if this value is true.}
  \item{resvec}{...}
}
\details{This routine tries to restore the temporal NDVI profile, i.e. separate true observations from noise. 
An algorithm creates a gaussian function fitted to NDVI composites.}
\value{The computed green-up date as Julian day (Day of Year).}
\author{Daniel Doktor, Maximilian Lange}

\examples{
values <- vector(mode="numeric", length=365)
values <- c(seq(from=1000, by=33, length=150), 
	rep(6000, 65), seq(from=6000, by=-33, length=150))
obsday <- 0
obsday <- gauss(values)
obsday
}

