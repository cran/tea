GH <-
function(data){
n=length(data)

Qn=function(k){

Tn=function(k){
  U=c();help=c();
  xstat = sort(data, decreasing = TRUE)
  for (i in 1:k) {
    U[i]=i*(log(xstat[i]) - log(xstat[i + 1]))
    help[i]=(k-2*i+1)*U[i]
  }
  Z=sum(help)
  N=mean(U)
  T=sqrt(3/k^3)*(Z/N)
  T
}

  start=k-floor(k/2)
  end=k+floor(k/2)
  y=start:end
  x=sapply(y,Tn)
  erg=1/(2*floor(k/2)+1)*sum(x^2)
  erg2=sqrt(erg)
  crit=1.25
  as.numeric(erg2>=crit)
}

kmax=floor(n/1.5)
i=1
while (Qn(i)==0){
  i=i+1
  if (i==kmax) break
}

u=sort(data,decreasing=TRUE)[i]

helphill=function (k) {
  xstat = sort(data, decreasing = TRUE)
  xihat = mean((log(xstat[1:k]) - log(xstat[k + 1])))
  xihat
}

ti=1/helphill(i)
list=list(k0=i,threshold=u,tail.index=ti)
list
}
