\name{analyzeBits}
\alias{analyzeBits}

\title{Bit Analyse}
\description{Analyses bits of a 16-bit integer}
\usage{analyzeBits(value, mode)}
\arguments{
  \item{value}{A 16-bit integer value}
  \item{mode}{If mode is one (default), the algorithm returns the most significant bit, 
if mode is two it returns the second significant bit and if mode is 3 the return value is the integer value of the last four bits}
}
\details{This routine analyses 16-bit integer values to get the indicators of MEDOKADS data}
\author{Daniel Doktor, Maximilian Lange}

\examples{
value <- -32768
res <- analyzeBits(value, mode=3)
res
}

