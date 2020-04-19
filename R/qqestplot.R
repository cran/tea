#' QQ-Estimator-Plot
#'
#' Plots the QQ-Estimator against the upper order statistics
#' @param data vector of sample data
#' @param kmin gives the minimal \code{k} for which the graph is plotted. Default ist set to \code{kmin=5}
#' @param conf.int \code{logical}. If FALSE (default) no confidence intervals are plotted
#' @details The QQ-Estimator is a Tail Index Estimator based on regression diagnostics. Assuming a Pareto tail behaviour of the data at hand a QQ-Plot of the theoretical quantiles of an exponential distribution against the empirical quantiles of the log-data should lead to a straight line above some unknown upper order statistic \code{k}. The slope of this line is an estimator for the tail index. Computing this estimator via linear regression for every \code{k} the plot should stabilize for the correct number of upper order statistics, denoted \code{k0} here.
#' @return The plot shows the values of the QQ-Estimator with respect to \code{k}. See references for more information.
#' @references Kratz, M. and Resnick, S.I. (1996). The QQ-estimator and heavy tails. \emph{Stochastic Models}, \bold{12(4)}, 699--724.
#' @examples
#' data(danish)
#' qqestplot(danish)
#' @export
qqestplot <-
function(data,kmin=5,conf.int=FALSE){
data=sort(data)
n=length(data)

QQplot=function(k){
  i=1:k
  xaxis=-log(1-(i/(k+1)))
  yaxis=log(data[n-k+i])
  SL=as.vector(lm(yaxis~xaxis)[[1]][2])
  return(1/SL)
}
x=1:n
y=sapply(x,QQplot)
plot(x[kmin:(n-1)],y[kmin:(n-1)],type="l",xlim=c(15,n),xlab="Order Statistics",ylab="qq-alpha")
if (conf.int == TRUE){
    confint.d = 1/y - (qnorm(0.975) * sqrt(2)/(y*sqrt(x)))
    confint.u = 1/y + (qnorm(0.975) * sqrt(2)/(y*sqrt(x)))
    lines(x[kmin:(n-1)],1/confint.d[kmin:(n-1)], lty = 2, col = "blue")
    lines(x[kmin:(n-1)],1/confint.u[kmin:(n-1)], lty = 2, col = "blue")
  }
}
