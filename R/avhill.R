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
