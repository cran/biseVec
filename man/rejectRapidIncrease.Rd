\name{rejectRapidIncrease}
\alias{rejectRapidIncrease}

\title{Rapid Increase Rejection}
\description{Rejects rapid increases of bise-algorithm resolutions}
\usage{rejectRapidIncrease(res)}
\arguments{
  \item{res}{NDVI values computed from the bise algorithm, values between 0 and 1 as a numeric vector.}
}
\value{A vector of numeric values between 0 and 1 without rapid increases (unnormal peaks).}
\author{Daniel Doktor, Maximilian Lange}

\examples{
values <- vector(mode="numeric", length=365)
values <- c(seq(from=0, by=0.01, length=50), 
	1, seq(from=0.5, by=-0.01, length=49))
res <- rejectRapidIncrease(values)
res
}

