% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qqestplot.R
\name{qqestplot}
\alias{qqestplot}
\title{QQ-Estimator-Plot}
\usage{
qqestplot(data, kmin = 5, conf.int = FALSE)
}
\arguments{
\item{data}{vector of sample data}

\item{kmin}{gives the minimal \code{k} for which the graph is plotted. Default ist set to \code{kmin=5}}

\item{conf.int}{\code{logical}. If FALSE (default) no confidence intervals are plotted}
}
\value{
The plot shows the values of the QQ-Estimator with respect to \code{k}. See references for more information.
}
\description{
Plots the QQ-Estimator against the upper order statistics
}
\details{
The QQ-Estimator is a Tail Index Estimator based on regression diagnostics. Assuming a Pareto tail behaviour of the data at hand a QQ-Plot of the theoretical quantiles of an exponential distribution against the empirical quantiles of the log-data should lead to a straight line above some unknown upper order statistic \code{k}. The slope of this line is an estimator for the tail index. Computing this estimator via linear regression for every \code{k} the plot should stabilize for the correct number of upper order statistics, denoted \code{k0} here.
}
\examples{
data(danish)
qqestplot(danish)
}
\references{
Kratz, M. and Resnick, S.I. (1996). The QQ-estimator and heavy tails. \emph{Stochastic Models}, \bold{12(4)}, 699--724.
}
