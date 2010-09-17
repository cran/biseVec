\name{createCoef}
\alias{createCoef}

\title{Savitzky-Golay Coefficient Creation}
\description{Creates coefficients for a Savitzky-Golay smoothing filter}
\usage{createCoef(window, degree)}
\arguments{
  \item{window}{The window size of the smoothing algorithm.}
  \item{degree}{Degree of the fitting polynomial. (2=quadratic filter, 4=quartic filter)}
}
\details{This routine computes the coefficients for a Savitzky-Golay filter 
based on the size of the moving window and the degree of the fitting polynomial}
\value{The computed coefficients for a Savitzky-Golay filter}
\references{Numerical Recipes: The Art of Scientific Computing, 1992, Chapter 14.8}
\author{Daniel Doktor, Maximilian Lange}

\examples{
window <- 7
degree <- 2
coef <- vector(mode="numeric", length=window)
coef <- createCoef(window, degree)
coef
}

