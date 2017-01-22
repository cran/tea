\name{avhill}
\alias{avhill}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Averaged Hill Plot
}
\description{
Plots an averaged version of the classical Hill Plot
}
\usage{
avhill(data, u = 2, kmin = 5, conf.int = FALSE)
}
\arguments{
  \item{data}{
vector of sample data
}
  \item{u}{
gives the amount of which the Hill estimator is averaged. Default ist set to \code{u=2}.
}
\item{kmin}{
gives the minimal \code{k} for which the graph is plotted. Default ist set to \code{kmin=5}.
}
  \item{conf.int}{
\code{logical}. If FALSE (default) no confidence intervals are plotted
}
}
\details{
The Averaged Hill Plot is a smoothed versionof the classical Hill Plot by taking the mean of values of the Hill estimator for subsequent \code{k}, i.e. upper order statistics. For more information see references. 
}
\value{
The normal black line gives the classical Hill Plot. The red dotted line is an averaged version that smoothes the Hill Plot by taking the mean of \code{k(u-1)} subsequent Hill estimations with respect to \code{k}. See references for more information.
}
\references{
Resnick, S. and Starica, C. (1997)
Smoothing the Hill estimator.
\emph{Advances in Applied Probability}, 271--293.
}
\examples{
data(danish)
avhill(danish)    
}
