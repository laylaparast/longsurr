\name{resam}
\alias{resam}
\title{
Resampling for Semiparametric Joint Linear Model
}
\description{
Resamples data for variance estimation for the semiparametric joint linear model estimator using weights
}
\usage{
resam(v, X, Time, Delta, obsT, Y)
}
\arguments{
  \item{v}{
resampling or perturbation weight, must be the same length of X
}
  \item{X}{
  numeric vector containing the treatment indicator for each observation with 1 for treated and 0 for control
}
  \item{Time}{
  numeric vector containing the observed event or censoring time for each observation 
}
  \item{Delta}{
  numeric vector containing the event indicator for each observation 

}
  \item{obsT}{
  numeric matrix containing the time that the surrogate marker was measured for each observation; number of rows should be equal to the number of observations and number of columns should be equal to the maximum number of surrogate markers measured. If the surrogate marker was not measured, the corresponding entry should be 0 or NA.

}
  \item{Y}{
numeric matrix containing the the surrogate marker measurements over time for each observation; number of rows should be equal to the number of observations and number of columns should be equal to the maximum number of surrogate markers measured. If the surrogate marker was not measured, as determined by the obsT entry, the Y at that time will be ignored.

}
}

\value{
Returns a numeric vector of resampled estimates.
}

\keyword{internal}
