\name{biseGauss}
\alias{biseGauss}

\title{Gaussian Function Creation trough Best Index Slope Extraction}
\description{Calculates the green-up day out of true NDVI values extracted from a gaussian function fitted to NDVI values as selected by a BISE-algorithm}
\usage{biseGauss(ndvi, slidperiod = 40, globalthreshold = 0, asym=FALSE, resvec=c(0.55,0.6,0.65,0.65,0.7))}
\arguments{
  \item{ndvi}{NDVI values computed from the daily surface reflectances, values between 0 and 10000 as a numeric vector.}
  \item{slidperiod}{Sliding Period of the BISE-algorithm, the default value is 40.}
  \item{globalthreshold}{The algorithm returns the green-up date calculated through a local threshold method if the value of mode is 0 or through a global threshold method (with a global threshold equal to the value of globalthreshold) if the value is greater than 0.}
  \item{asym}{An asymmetric Gaussian function will be created if this value is true.}
  \item{resvec}{...}
}
\details{This routine tries to restore the temporal NDVI profile, i.e. separate true observations from noise. 
An algorithm creates a gaussian function fitted to NDVI values as selected by a BISE-algorithm (details in Viovy et. al (1992)).
The sliding period default is 40 days, optimised for the area of Central Europe. 
This might has to be modified when vegetation dynamics of other climates/habitats are analysed.}
\value{The computed green-up date as Julian day (Day of Year).}
\references{Viovy et. al (1992). The Best Index Slope Extraction (BISE) - a method for reducing noise in NDVI time-series. International Journal of Remote Sensing, 13, 1585-1590.}
\author{Daniel Doktor, Maximilian Lange}

\examples{
values <- vector(mode="numeric", length=365)
values <- c(seq(from=1000, by=33, length=150), 
	rep(6000, 65), seq(from=6000, by=-33, length=150))
obsday <- 0
obsday <- biseGauss(values)
obsday
}

