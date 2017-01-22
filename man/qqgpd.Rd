\name{qqgpd}
\alias{qqgpd}
\title{
QQ-Plot against the generalized Pareto distribution for given number of exceedances
}
\description{
Plots the empirical observations above a given threshold against the theoretical quantiles of a generalized Pareto distribution.
}
\usage{
qqgpd(data, nextremes, scale, shape)
}
\arguments{
  \item{data}{vector of sample data}
  \item{nextremes}{number of exceedances}
  \item{scale}{scale parameter of the GPD}
  \item{shape}{shape parameter of the GPD}
}
\details{
If the fitted GPD model provides a reasonable approximation of the underlying sample data the empirical and theoretical quantiles should coincide. So plotting them against each other should result in a straight line. Deviations from that line speak for a bad model fit and against a GPD assumption. 
}
\value{
The straight red line gives the line of agreement. The dashed lines are simulated 95 percent confidence intervals. Therefor the fitted GPD model is simulated 1000 times using Monte Carlo. The sample size of each simulation equals the number of exceedances. 
}
\examples{
data=rexp(1000) #GPD with scale=1, shape=0
qqgpd(data,1000,1,0)
}
