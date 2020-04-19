#' Averaged Hill Plot
#'
#' Plots an averaged version of the classical Hill Plot
#' @param data vector of sample data
#' @param u gives the amount of which the Hill estimator is averaged. Default ist set to \code{u=2}.
#' @param kmin gives the minimal \code{k} for which the graph is plotted. Default ist set to \code{kmin=5}.
#' @param conf.int \code{logical}. If FALSE (default) no confidence intervals are plotted
#' @details The Averaged Hill Plot is a smoothed version of the classical Hill Plot by taking the mean of values of the Hill estimator for subsequent \code{k}, i.e. upper order statistics. For more information see references. 
#' @return The normal black line gives the classical Hill Plot. The red dotted line is an averaged version that smoothes the Hill Plot by taking the mean of \code{k(u-1)} subsequent Hill estimations with respect to \code{k}. See references for more information.
#' @references Resnick, S. and Starica, C. (1997). Smoothing the Hill estimator. \emph{Advances in Applied Probability}, 271--293.
#' @examples
#' data(danish)
#' avhill(danish) 
#' @export
avhill <-
function(data,u=2,kmin=5, conf.int=FALSE){
n=length(data)

i=1:(n-1)
xstat=sort(data,decreasing=TRUE)
h=(cumsum(log(xstat[i]))/i)-log(xstat[i+1]) 

plot(i[kmin:(n-1)],1/h[kmin:(n-1)],type="l",xlab="Order Statistics",ylab="Tail Index") 
confint.d=h-(qnorm(0.975)*h/sqrt(i))
confint.u=h+(qnorm(0.975)*h/sqrt(i))

x=1:floor(n/u)
y=c();
for (k in 1:length(x)) { 
  y[k]=mean(h[(k+1):(u*k)])
}
lines(x,1/y,type="l",col="red")
scal=sqrt( 2/(u-1)*( 1-(log(u)/(u-1)) )  )
confint2.d=y-(qnorm(0.975)*y*scal/sqrt(x))
confint2.u=y+(qnorm(0.975)*y*scal/sqrt(x))

if (conf.int==TRUE){
  lines(i[kmin:(n-1)],1/confint.d[kmin:(n-1)],lty=2,col="blue")
  lines(i[kmin:(n-1)],1/confint.u[kmin:(n-1)],lty=2,col="blue")
  lines(x[kmin:(n-1)],1/confint2.d[kmin:(n-1)],lty=2,col="green")
  lines(x[kmin:(n-1)],1/confint2.u[kmin:(n-1)],lty=2,col="green")
 }
}
