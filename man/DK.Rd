\name{DK}
\alias{DK}
\title{
A Bias-based procedure for Choosing the Optimal Sample Fraction
}
\description{
An Implementation of the procedure proposed in Drees & Kaufmann (1998) for selecting the optimal sample fraction in tail index estimation.
}
\usage{
DK(data, r=1)
}
\arguments{
  \item{data}{vector of sample data}
  \item{r}{tuning parameter for the stopping criterion. \code{default} is set to \code{1}. Change only if recommended by the output.}
}
\details{
The procedure proposed in Drees & Kaufmann (1998) is based on bias reduction. A stopping criterion with respect to \code{k} is implemented to find the optimal tail fraction, i.e. \code{k/n} with \code{k} the optimal number of upper order statistics. This number, denoted \code{k0} here, is equivalent to the number of extreme values or, if you wish, the number of exceedances in the context of a POT-model like the generalized Pareto distribution. \code{k0} can then be associated with the unknown threshold \code{u} of the GPD by choosing \code{u} as the \code{n-k0}th upper order statistic. If the above mentioned stopping criterion exceedes a certain value \code{r}, the bias of the assumed extreme model has become prominent and therefore \code{k} should not be chosen higher. For more information see references.
}
\value{  \item{second.order.par }{gives an estimation of the second order parameter \code{rho}.}
  \item{k0}{optimal number of upper order statistics, i.e. number of exceedances or data in the tail}
  \item{threshold}{the corresponding threshold}
  \item{tail.index}{the corresponding tail}
}
\references{
Drees, H. and Kaufmann, E. (1998)
Selecting the optimal sample fraction in univariate extreme value estimation.
\emph{Stochastic Processes and their Applications}, \bold{75(2)}, 149--172.
}
\examples{
data(danish)
DK(danish)
}