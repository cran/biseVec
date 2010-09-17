\name{SavGol}
\alias{SavGol}

\title{Savitzky-Golay Function Creation}
\description{Calculates the green-up day out of true NDVI values extracted from a Savitzky-Golay function fitted to NDVI composites}
\usage{SavGol(ndvi, window=7, degree=2, smoothing=10, globalthreshold = 0, resvec=c(0.55,0.6,0.65,0.65,0.7))}
\arguments{
  \item{ndvi}{NDVI values (NDVI composites) computed from the daily surface reflectances, values between 0 and 10000 as a numeric vector.}
  \item{window}{The window size of the smoothing algorithm.}
  \item{degree}{Degree of the fitting polynomial. (2=quadratic filter, 4=quartic filter)}
  \item{smoothing}{How often the values will be smoothed.}
  \item{globalthreshold}{The algorithm returns the green-up date calculated through a local threshold method if the value of mode is 0 or through a global threshold method (with a global threshold equal to the value of globalthreshold) if the value is greater than 0.}
  \item{resvec}{...}
}
\details{This routine tries to restore the temporal NDVI profile, i.e. separate true observations from noise. 
An algorithm creates a Savitzky-Golay function fitted to NDVI composites.}
\value{The computed green-up date as Julian day (Day of Year).}
\references{Numerical Recipes: The Art of Scientific Computing, 1992, Chapter 14.8}
\author{Daniel Doktor, Maximilian Lange}

\examples{
values <- vector(mode="numeric", length=365)
values <- c(seq(from=1000, by=33, length=150), 
	rep(6000, 65), seq(from=6000, by=-33, length=150))
obsday <- 0
obsday <- SavGol(values)
obsday
}

