RT <-
function(data,beta=0,kmin=2){
n=length(data)

x1=sort(data, decreasing=TRUE)
i=1:(n-1)
j=i^beta
h=(cumsum(log(x1[i]))/i)-log(x1[i+1])

med=c();
for (i in 1:(n-1)){
  med[i]=median(h[1:i])
}
erg1=c(); erg2=c()
for (k in 1:(n-1)){
erg1[k]=sum(j[1:k]*abs(h[1:k]-med[k]))/k
erg2[k]=sum(j[1:k]*(h[1:k]-h[k])^2)/k
}
rt1=which.min(erg1[kmin:(n-1)])
rt2=which.min(erg2[kmin:(n-1)])
list=list(k0=c(rt1,rt2),threshold=c(x1[rt1],x1[rt2]),tail.index=c(1/h[rt1],1/h[rt2]))
list
}
