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
