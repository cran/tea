\name{dAMSE}
\alias{dAMSE}
\title{
Minimizing the AMSE of the Hill estimator with respect to k
}
\description{
Gives the optimal number of upper order statistics \code{k} for the Hill estimator by minimizing the AMSE-criterion.
}
\usage{
dAMSE(data)
}
\arguments{
  \item{data}{vector of sample data}
}
\details{
The optimal number of upper order statistics is equivalent to the number of extreme values or, if you wish, the number of exceedances in the context of a POT-model like the generalized Pareto distribution. This number is identified by minimizing the AMSE criterion with respect to \code{k}. The optimal number, denoted \code{k0} here, can then be associated with the unknown threshold \code{u} of the GPD by choosing \code{u} as the \code{n-k0}th upper order statistic. For more information see references.
}
\value{
  \item{second.order.par}{gives an estimation of the second order parameter \code{beta} and \code{rho}.}
  \item{k0}{optimal number of upper order statistics, i.e. number of exceedances or data in the tail}
  \item{threshold}{the corresponding threshold}
  \item{tail.index}{the corresponding tail index}
}
\references{
Caeiro, J. and Gomes, M.I. (2016)
Threshold selection in extreme value analysis.
\emph{Extreme Value Modeling and Risk Analysis:Methids and Applications}, 69--86.
}
\examples{
data(danish)
dAMSE(danish)
}
