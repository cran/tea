althill <-
function(data,u=2,kmin=5,conf.int=FALSE){
n=length(data)

hill=function (k1){ 
  xstat = sort(data, decreasing = TRUE)
  xihat = mean(log(xstat[1:k1])) - log(xstat[k1 + 1])
  xihat
}

avhill=function(k){ 
  frac=1/((u-1)*k)
  H_pn=sapply((k+1):(u*k),hill)
  avHill=frac*sum(H_pn)
  avHill
}

theta=seq(0,1,0.01)
help=ceiling(n^theta)
alt1=sapply(help,hill) 
plot(theta[kmin:(n-1)],1/alt1[kmin:(n-1)],type="l",xlab="theta",ylab="Tail Index") 
confint.d=alt1-(qnorm(0.975)*alt1/sqrt(help))
confint.u=alt1+(qnorm(0.975)*alt1/sqrt(help))


alt2=sapply(help,avhill) 
lines(theta[kmin:(n-1)],1/alt2[kmin:(n-1)],col="red") 
scal=sqrt( 2/(u-1)*( 1-(log(u)/(u-1)) )  )
confint2.d=alt2-(qnorm(0.975)*alt2*scal/sqrt(help))
confint2.u=alt2+(qnorm(0.975)*alt2*scal/sqrt(help))

if(conf.int==TRUE){
  lines(theta[kmin:(n-1)],1/confint.d[kmin:(n-1)],lty=2,col="blue")
  lines(theta[kmin:(n-1)],1/confint.u[kmin:(n-1)],lty=2,col="blue")
  lines(theta[kmin:(n-1)],1/confint2.d[kmin:(n-1)],lty=2,col="green")
  lines(theta[kmin:(n-1)],1/confint2.u[kmin:(n-1)],lty=2,col="green")
}
}
