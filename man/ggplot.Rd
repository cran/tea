\name{ggplot}
\alias{ggplot}
\title{
Gerstengarbe Plot
}
\description{
Performs a sequential Mann-Kendall Plot also known as Gerstengarbe Plot.
}
\usage{
ggplot(data, nexceed = min(data) - 1)
}
\arguments{
  \item{data}{vector of sample data}
  \item{nexceed}{number of exceedances. Default is the minimum of the data to make sure the whole dataset is considered.}
}
\details{
The Gerstengarbe Plot, referring to Gerstengarbe and Werner (1989), is a sequential version of the Mann-Kendall-Test. This test searches for change points within a time series. This method is adopted for finding a threshold in a POT-model. The basic idea is that the differences of order statistics of a given dataset behave different between the body and the tail of a heavy-tailed distribution. So there should be a change point if the POT-model holds. 
To identify this change point the sequential test is done twice, for the differences from start to the end of the dataset and vice versa. The intersection point of these two series can then be associated with the change point of the sample data. For more informations see references.
}
\value{
  \item{k0}{optimal number of upper order statistics, i.e. the change point of the dataset}
  \item{threshold}{the corresponding threshold}
  \item{tail.index}{the corresponding tail index}
}
\section{Authors}{
Ana Cebrian

Johannes Ossberger
}
\section{Acknowledgements}{
Great thanks to A. Cebrian for providing a basic version of this code. 
}
\references{
Gerstengarbe, F.W. and Werner, P.C. (1989)
A method for statistical definition of extreme-value regions and their application to meteorological time series.
\emph{Zeitschrift fuer Meteorologie}, \bold{39(4)}, 224--226.

Cebrian, A., and Denuit, M. and Lambert, P. (2003)
Generalized pareto fit to the society of actuaries large claims database.
\emph{North American Actuarial Journal}, \bold{7(3)}, 18--36.
} 
\examples{
data(danish)
ggplot(danish)
}
