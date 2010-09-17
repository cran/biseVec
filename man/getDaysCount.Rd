\name{getDaysCount}
\alias{getDaysCount}

\title{Number of Days}
\description{Returns the number of days of the year}
\usage{getDaysCount(date)}
\arguments{
	\item{date}{Date YYMMDD as integer}
}
\value{The number of days of the year in the date}
\author{Daniel Doktor, Maximilian Lange}

\examples{
	date <- 891208
	days <- 0
	days <- getDaysCount(date)
	days
}