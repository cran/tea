\name{althill}
\alias{althill}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Alternative Hill Plot
}
\description{
Plots the Alternative Hill Plot and an averaged version of it against the upper order statistics.
}
\usage{
althill(data, u = 2, kmin = 5, conf.int = FALSE)
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
The Alternative Hill Plot is just a normal Hill Plot scaled to the \code{[0,1]} interval which can make interpretation much easier. See references for more information.
}
\value{
The normal black line gives a simple Hill Plot scaled to \code{[0,1]}. The red dotted line is an averaged version that smoothes the Hill Plot by taking the mean of \code{k(u-1)} subsequent Hill estimations with respect to \code{k}. See references for more information.
}
\references{
Resnick, S. and Starica, C. (1997)
Smoothing the Hill estimator.
\emph{Advances in Applied Probability}, 271--293.
}
\examples{
data=rexp(500)
althill(data)  
}
