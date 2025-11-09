\name{data_sjm}
\alias{data_sjm}
\docType{data}
\title{
Example data for semiparametric joint estimation functions
}
\description{
Simulated example data for semiparametric joint Estimation functions
}
\usage{data("data_sjm")}
\format{

 A list with 200 observations on the following:
  \describe{
    \item{\code{delta}}{numeric vector containing the event indicator for each observation}
    \item{\code{obsT}}{numeric matrix containing the time that the surrogate marker was measured for each observation; number of rows is equal to the number of observations (200) and number of columns is equal to the maximum number of surrogate markers measured (15)}
      \item{\code{Y}}{numeric matrix containing the surrogate marker measurements over time for each observation; same dimension as obsT}
     \item{\code{Time}}{numeric vector containing the observed event or censoring time for each observation }
       \item{\code{Treatment}}{numeric vector containing the treatment indicator for each observation with 1 for treated and 0 for control}
  }
}

\examples{
data(data_sjm)
names(data_sjm)
}
\keyword{datasets}
