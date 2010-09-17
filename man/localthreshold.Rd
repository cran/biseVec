\name{localThresGreenup}
\alias{localThresGreenup}

\title{Green-Up Day Calculation by Local Threshold Method}
\description{Calculates the green-up day out of fitted NDVI values by a local threshold method}
\usage{localThresGreenup(model, localthreshold)}
\arguments{
  \item{model}{The model for which the green-up day shut be calculated.}
  \item{localthreshold}{This value determines the percentage between minimum and maximum of model at which the green-up day will be set. }
}
\details{This routine calculates the green-up day out of a model of fitted NDVI values by a local threshold method.}
\value{The computed green-up date as Julian day (Day of Year).}
\author{Daniel Doktor, Maximilian Lange}

\examples{\dontrun{
#model with fitted NDVI values
model

greenup <- localThresGreenup(model, 0.6)
}
}

