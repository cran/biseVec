\name{globalThresGreenup}
\alias{globalThresGreenup}

\title{Green-Up Day by Global Threshold}
\description{Calculates the green-up day out of fitted NDVI values by a global threshold method}
\usage{globalThresGreenup(model, globalthreshold)}
\arguments{
  \item{model}{The model for which the green-up day shut be calculated.}
  \item{globalthreshold}{This value determines the threshold of this method. If the model reaches this threshold, the day will returned as the green-up day. }
}
\details{This routine calculates the green-up day out of a model of fitted NDVI values by a global threshold method.}
\value{The computed green-up date as Julian day (Day of Year).}
\author{Daniel Doktor, Maximilian Lange}

\examples{\dontrun{
#model with fitted NDVI values
model

greenup <- globalThresGreenup(model, 0.6)
}
}

