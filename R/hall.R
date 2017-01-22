hall <-
function(data,B=1000,epsilon=0.955,kaux=2*sqrt(length(data))){
n=length(data)
n1=floor(n^epsilon)

helphill=function (k) {
  xstat = sort(data, decreasing = TRUE)
  xihat = mean((log(xstat[1:k]) - log(xstat[k + 1])))
  xihat
}

help=helphill(kaux)

mse=matrix(nrow=B,ncol=n1-1)
for (l in 1:B){
  i=1:(n1-1)
  x1=sample(data,n1,replace=TRUE)
  x1=sort(x1,decreasing=TRUE)
  h=(cumsum(log(x1[i]))/i)-log(x1[i+1]) 
  mse[l,]=h-help
}
mse=mse^2
msestar=colMeans(mse)
k1star=which.min(msestar)

k0star=floor(k1star*(n/n1)^(2/3))
u=sort(data,decreasing=TRUE)[k0star]
ti=1/helphill(k0star)
list=list(k0=k0star,threshold=u,tail.index=ti)
list
}
