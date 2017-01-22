HW <-
function(data){
n=length(data)
sigma=0.5
tau1=0.9
tau2=0.95
s=floor(n^sigma)
t1=floor(n^tau1)
t2=floor(n^tau2)

helphill=function (k) {
  xstat = sort(data, decreasing = TRUE)
  xihat = mean((log(xstat[1:k]) - log(xstat[k + 1])))
  xihat
}

Z=(1/helphill(t1))-(1/helphill(s))
N=(1/helphill(t2))-(1/helphill(s))
ZN=abs(Z/N)
NN=t1/t2
rho=abs(log(ZN)/log(NN)) 

frac=(1/helphill(t1)-1/helphill(s))/helphill(s)
fac=(1/sqrt(2*rho))*((n/t1)^rho)
Exp=2/(2*rho+1)
lambda=(abs(fac*frac))^Exp 

k0star=floor(lambda*n^((2*rho)/(2*rho+1)))
u=sort(data,decreasing=TRUE)[k0star]
ti=1/helphill(k0star)
list=list(sec.order.par=-rho,k0=k0star,threshold=u,tail.index=ti)

c1=-sigma/(2*(1-tau1))
c2=-sigma/(2*(1-sigma))
if (rho < c1 && rho > c2) print("Warning: check consistency! Use different tuning parameters")
list
}
