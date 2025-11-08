\name{presmooth_data}
\alias{presmooth_data}
\title{Pre-smooth sparse longitudinal data}
\usage{
presmooth_data(obs_data, ...)
}
\arguments{
\item{obs_data}{data.frame or tibble containing the observed data, with columns \code{id} identifying the individual measured, \code{tt} identifying the time of the observation, \code{x} the value of the surrogate at time \code{tt}, and \code{a} indicating 1 for treatment arm and 0 for control arm.}
\item{...}{additional arguments passed on to \code{fpca}}
}
\value{
list containing matrices \code{X_t} and \code{X_c}, which are the smoothed surrogate values for the treated and control groups, respectively, for use in downstream analyses
}
\description{
Pre-smooth sparse longitudinal data
}
\examples{
library(dplyr)
data(full_data)
obs_ds <- group_by(full_data, id) 
obs_data <- sample_n(obs_ds, 5)
obs_data <- ungroup(obs_data)

head(obs_data)
presmooth_X <- presmooth_data(obs_data)
}