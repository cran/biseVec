\name{biseDSigC}
\alias{biseDSigC}

\title{Double Sigmoid Function Creation trough Best Index Slope Extraction with C-Subroutine}
\description{Calculates the green-up day out of true NDVI values extracted from a double sigmoid function fitted to NDVI values (with C-Subroutine instead of nls()) as selected by a BISE-algorithm}
\usage{biseDSigC(ndvi, slidperiod = 40, globalthreshold = 0, localthreshold = 0.55, draw=FALSE, color="blue",resvec=c(0.55,0.6,0.65,0.65,0.7))}
\arguments{
  \item{ndvi}{NDVI values computed from the daily surface reflectances, values between 0 and 10000 as a numeric vector.}
  \item{slidperiod}{Sliding Period of the BISE-algorithm, the default value is 40.}
  \item{globalthreshold}{The algorithm sets the green-up day at the inflection point of the function if the value is less than 0. 
If the value is 0, the green-up day is calculated through a local threshold method and for a value greater than 0 through a global threshold method (with a global threshold equal to the value of globalthreshold). }
  \item{localthreshold}{If globalthreshold is 0, this value determines the percentage between minimum and maximum of the function at which the green-up day will be set. The default value is 0.55.}
  \item{draw}{A boolean value which determines if the function should be drawn into an open plot.}
  \item{color}{If draw is true, the function is drawn with this color.}
  \item{resvec}{...}
}
\details{This routine tries to restore the temporal NDVI profile, i.e. separate true observations from noise. 
An algorithm creates a double sigmoid function fitted to NDVI values as selected by a BISE-algorithm (details in Viovy et. al (1992)).
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
obsday <- biseDSigC(values)
obsday
}

