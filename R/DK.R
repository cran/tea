DK <-
function(data,r=1){
n=length(data)

helphill=function (k) {
  xstat = sort(data, decreasing = TRUE)
  xihat = mean((log(xstat[1:k]) - log(xstat[k + 1])))
  xihat
}

gamma.tilde=helphill(2*sqrt(n))
r.n=r*2.5*gamma.tilde*n^0.25

kn=function(k){
  i=1:k
  x=sort(data,decreasing=TRUE)
  h=(cumsum(log(x[i]))/i)-log(x[i+1]) 
  h2=sqrt(i)*abs(h-h[k])
  Max=max(which(h2==max(h2)))
  Max
}
khelp=sapply(1:(n-1),kn)
khelp2=khelp>r.n
kbar_r.n=min(which(khelp2==1))

if (kbar_r.n==Inf) {
  print("Warning: no k_n>r_n found. Use smaller r_n. Try r=0.9")
} else {
  xi=0.7
  lambda=0.6
  r.n.xi=r.n^xi
  
  khelp3=khelp>r.n.xi
  kbar_r.n.xi=min(which(khelp3==1))
  
  Z=log(kn(floor(lambda*kbar_r.n.xi)))-log(kn(kbar_r.n.xi))
  rho=-(Z/log(lambda))+0.5
    
  if (rho>=0) {
    print("Warning: rho>0. Method fails here.")
  } else {
    Exp=1/(1-2*rho)
    kk=(kbar_r.n.xi/(kbar_r.n^xi))^(1/(1-xi))
    k0star=(((1-2*rho)^(1/rho))*((-2*rho*gamma.tilde^2)^Exp)*kk)
    k0star=floor(k0star)
    ti=1/helphill(k0star)
    u=sort(data,decreasing=TRUE)[k0star]
    list=list(sec.order.par=rho,k0=k0star,threshold=u,tail.index=ti)
    list
  }
}
}
