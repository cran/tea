danielsson <-
function(data,B=500,epsilon=0.9){
n=length(data)
n1=floor(n^epsilon)
n2=(n1^2)/n

Qn1=function (k) { 
  xstat = sort(x1, decreasing = TRUE)
  xihat = mean((log(xstat[1:k]) - log(xstat[k + 1]))^2)
  xihat2 = mean((log(xstat[1:k]) - log(xstat[k + 1])))
  xihat-(2*xihat2^2)
}
Qn2=function (k) { 
  xstat = sort(x2, decreasing = TRUE)
  xihat = mean((log(xstat[1:k]) - log(xstat[k + 1]))^2)
  xihat2 = mean((log(xstat[1:k]) - log(xstat[k + 1])))
  xihat-(2*xihat2^2)
}

qn1=matrix(nrow=B,ncol=n1-1)
qn2=matrix(nrow=B,ncol=n2-1)
for (l in 1:B){
  x1=sample(data,n1,replace=TRUE)
  x2=sample(data,n2,replace=TRUE)
  qn1[l,]=sapply(1:(n1-1),Qn1)
  qn2[l,]=sapply(1:(n2-1),Qn2)
}
qn1=qn1^2
qn2=qn2^2
qn1star=colMeans(qn1)
qn2star=colMeans(qn2)
k1star=which.min(qn1star)
k2star=which.min(qn2star)

Exp=(log(n1)-log(k1star))/(log(n1))
Z=(log(k1star))^2
N=(2*log(n1)-log(k1star))^2
k0star=floor(k1star^2/k2star*((Z/N)^Exp))+1
u=sort(data,decreasing=TRUE)[k0star]
rho=log(k1star)/(-2*log(n1)+2*log(k1star))

helphill=function (k) {
  xstat = sort(data, decreasing = TRUE)
  xihat = mean((log(xstat[1:k]) - log(xstat[k + 1])))
  xihat
}

ti=1/helphill(k0star)
list=list(sec.order.par=rho,k0=k0star,threshold=u,tail.index=ti)
list
}
