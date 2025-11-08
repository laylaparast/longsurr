\name{estimate_surrogate_value}
\alias{estimate_surrogate_value}
\title{Estimate the surrogate value of a longitudinal marker}
\usage{
estimate_surrogate_value(y_t, y_c, X_t, X_c, method = c("gam", "linear",
  "kernel"), k = 3, var = FALSE, bootstrap_samples = 50, alpha = 0.05)
}
\arguments{
\item{y_t}{vector of n1 outcome measurements for treatment group}

\item{y_c}{vector of n0 outcome measurements for control or reference group}

\item{X_t}{n1 x T matrix of longitudinal surrogate measurements for treatment group, where T is the number of time points}

\item{X_c}{n0 x T matrix of longitudinal surrogate measurements for control or reference group, where T is the number of time points}

\item{method}{method for dimension-reduction of longitudinal surrogate, either 'gam', 'linear', or 'kernel'}

\item{k}{number of eigenfunctions to use in semimetric}

\item{var}{logical, if TRUE then standard error estimates and confidence intervals are provided}

\item{bootstrap_samples}{number of bootstrap samples to use for standard error estimation, used if var = TRUE, default is 50}

\item{alpha}{alpha level, default is 0.05}
}
\value{
a tibble containing estimates of the treatment effect (Deltahat), the residual treatment effect (Deltahat_S), and the proportion of treatment effect explained (R); if var = TRUE, then standard errors of Deltahat_S and R are also provided (Deltahat_S_se and R_se), and quantile-based 95\% confidence intervals for Deltahat_S and R are provided (Deltahat_S_ci_l [lower], Deltahat_S_ci_h [upper], R_ci_l [lower], R_ci_u [upper])
}
\description{
Estimate the surrogate value of a longitudinal marker
}
\references{
Agniel D and Parast L (2021). Evaluation of Longitudinal Surrogate Markers. Biometrics, 77(2): 477-489.
}
\examples{
library(dplyr)
data(full_data)


wide_ds <- full_data \%>\% 
dplyr::select(id, a, tt, x, y) \%>\%
tidyr::spread(tt, x) 

wide_ds_0 <- wide_ds \%>\% filter(a == 0)
wide_ds_1 <- wide_ds \%>\% filter(a == 1)
X_t <- wide_ds_1 \%>\% dplyr::select(`-1`:`1`) \%>\% as.matrix
y_t <- wide_ds_1 \%>\% pull(y)
X_c <- wide_ds_0 \%>\% dplyr::select(`-1`:`1`) \%>\% as.matrix
y_c <- wide_ds_0 \%>\% pull(y)

estimate_surrogate_value(y_t = y_t, y_c = y_c, X_t = X_t, X_c = X_c, 
method = 'gam', var = FALSE)
estimate_surrogate_value(y_t = y_t, y_c = y_c, X_t = X_t, X_c = X_c, 
method = 'linear', var = TRUE, bootstrap_sample = 50)
}
