\name{sjm_nl_estimate}
\alias{sjm_nl_estimate}
\title{
Semiparametric Joint Modeling of the Treatment Effect on a Longitudinal Surrogate with a Nonlinear Model}
\description{
Semiparametric joint modeling of the treatment effect on a longitudinal surrogate using both a Cox proportional hazards model and a splines-based model
}
\usage{
sjm_nl_estimate(X, Time, Delta, obsT, Y, gap_time = 0.1, n.resample = 100, var = FALSE)
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
\item{gap_time}{
number indicating gap time for slope estimation; default is 0.1.
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
  \item{est}{estimated hazard ratio from the Cox model}
  \item{est_t}{vector of estimated treatment effect on the slope of the surrogate marker i.e., the surrogate marker trajectory, on a grid constructed from the given gap time}
  \item{t_grid}{vector of grid times corresponding to the returned estimates}
  \item{SE_est}{if var is TRUE, standard error estimate of the hazard ratio}
   \item{SE_est_t}{if var is TRUE, standard error estimate of the estimated treatment effect on the slope of the surrogate marker}
  \item{CI_lower_est}{if var is TRUE, lower bound of the 95\% confidence interval for the hazard ratio}
  \item{CI_upper_est}{if var is TRUE, upper bound of the 95\% confidence interval for the hazard ratio}
    \item{CI_lower_est_t}{if var is TRUE, lower bound of the 95\% confidence interval for the treatment effect on the slope of the surrogate marker}
  \item{CI_upper_est_t}{if var is TRUE, upper bound of the 95\% confidence interval for the treatment effect on the slope of the surrogate marker}

}
\references{
Wang X, Zhou J, Parast L, Greene T (2025). Semiparametric Joint Modeling to Estimate
the Treatment Eﬀect on a Longitudinal Surrogate with Application to Chronic Kidney
Disease Trials. Biometrics, 81(3): ujaf104.
}
\author{
Xuan Wang
}

\examples{
data(data_sjm)
sjm_nl_estimate(X=data_sjm$Treatment, Time = data_sjm$Time, 
Delta = data_sjm$delta, obsT = data_sjm$obsT, Y = data_sjm$Y, gap_time=0.2)

\donttest{sjm_nl_estimate(X=data_sjm$Treatment, Time = 
data_sjm$Time, Delta = data_sjm$delta, obsT = data_sjm$obsT, 
Y = data_sjm$Y, gap_time = 0.2, n.resample=5, var=TRUE)}

}