\name{sjm_linear_estimate}
\alias{sjm_linear_estimate}
\title{
Semiparametric Joint Modeling of the Treatment Effect on a Longitudinal Surrogate with a Linear Model}
\description{
Semiparametric joint modeling of the treatment effect on a longitudinal surrogate using both a Cox proportional hazards model and linear model
}
\usage{
sjm_linear_estimate(X, Time, Delta, obsT, Y, n.resample=100, var = FALSE)
}
\arguments{
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
  \item{n.resample}{
number of resampled estimates used for variance estimation; default is 100.
}
  \item{var}{
logical indicating whether the user would like variance estimates and confidence intervals; default is FALSE.
}
}

\value{
A list of estimates is returned:
  \item{est}{vector of point estimates where the first entry is the hazard ratio from the Cox model, the second entry is the estimated treatment effect on the surrogate marker at baseline, and the third entry is the estimated treatment on the slope of the surrogate marker i.e., the surrogate marker trajectory}
  \item{SE}{if var is TRUE, a vector of standard error estimates corresponding to the returned point estimates}
  \item{CI_lower}{if var is TRUE, a vector of estimates for the lower bound of the 95\% confidence interval for the quantities corresponding to the returned point estimates}
  \item{CI_upper}{if var is TRUE, a vector of estimates for the upper bound of the 95\% confidence interval for the quantities corresponding to the returned point estimates}
}
\references{
Wang X, Zhou J, Parast L, Greene T (2025). Semiparametric Joint Modeling to Estimate
the Treatment Effect on a Longitudinal Surrogate with Application to Chronic Kidney
Disease Trials. Biometrics, 81(3): ujaf104.
 
}
\author{
Xuan Wang}

\examples{
data(data_sjm)
sjm_linear_estimate(X=data_sjm$Treatment, Time = data_sjm$Time, 
Delta = data_sjm$delta, obsT = data_sjm$obsT, Y = data_sjm$Y)

\donttest{sjm_linear_estimate(X=data_sjm$Treatment, Time = 
data_sjm$Time, Delta = data_sjm$delta, obsT = data_sjm$obsT, 
Y = data_sjm$Y, n.resample=5, var=TRUE)}

}