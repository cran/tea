Himp <-
function(data,B=1000,epsilon=0.955){
n=length(data)
n1=floor(n^epsilon)

helphill=function (k,j) {
  xstat = sort(data, decreasing = TRUE)
  xihat = mean((log(xstat[1:k]) - log(xstat[k + 1]))^j)
  xihat
}

DAMSE=function(data){
  n=length(data)
  k=c(floor(n^0.995),floor(n^0.999))
  
  M11=helphill(k[1],1)
  M12=helphill(k[1],2)
  M13=helphill(k[1],3)
  
  M21=helphill(k[2],1)
  M22=helphill(k[2],2)
  M23=helphill(k[2],3) 
  
  #tau=1
  W.k1.t1=(M11-sqrt(M12/2))/(sqrt(M12/2)-(M13/6)^(1/3))
  W.k2.t1=(M21-sqrt(M22/2))/(sqrt(M22/2)-(M23/6)^(1/3))
  rho.k1.t1=-abs(3*(W.k1.t1-1)/(W.k1.t1-3))
  rho.k2.t1=-abs(3*(W.k2.t1-1)/(W.k2.t1-3))
  
  #tau=0
  W.k1.t0=(log(M11)-0.5*log(M12/2))/(0.5*log(M12/2)-log(M13/6)/3)
  W.k2.t0=(log(M21)-0.5*log(M22/2))/(0.5*log(M22/2)-log(M23/6)/3)
  rho.k1.t0=-abs(3*(W.k1.t0-1)/(W.k1.t0-3))
  rho.k2.t0=-abs(3*(W.k2.t0-1)/(W.k2.t0-3))
  
  chi.t1=median(c(rho.k1.t1,rho.k2.t1))
  chi.t0=median(c(rho.k1.t0,rho.k2.t0))
  I.t1=c((rho.k1.t1-chi.t1)^2,(rho.k2.t1-chi.t1)^2)
  I.t0=c((rho.k1.t0-chi.t0)^2,(rho.k2.t0-chi.t0)^2)
  
  if (sum(I.t0)<=sum(I.t1)) { 
    tau=0; rho=rho.k2.t0
  } else {
    tau=1; rho=rho.k2.t1
  }
  
  U=c();
  xstat = sort(data, decreasing = TRUE)
  for (i in 1:k[2]) {
    U[i]=i*(log(xstat[i]) - log(xstat[i + 1]))
  }
  
  i=1:k[2]
  dk=mean((i/k[2])^(-rho))
  
  Dk=function(a){
    D=mean((i/k[2])^(-a)*U[i])
  }
  
  beta=(k[2]/n)^rho*(dk*Dk(0)-Dk(rho))/(dk*Dk(rho)-Dk(2*rho)) 
  erg=c(beta,rho)
  erg
}
secorder=DAMSE(data)

khelp1=2:(n1-1)
tn1=matrix(nrow=B,ncol=n1-2)
for (l in 1:B){
  x1=sample(data,n1,replace=TRUE)
  x1=sort(x1,decreasing=TRUE)
  i1=2:(n1-1)
  m1=floor(i1/2)
  Tkn1.1=(cumsum(log(x1[i1]))/i1)-log(x1[i1+1]) 
  Tkn1.2=(cumsum(log(x1[m1]))/i1)-log(x1[m1+1])
  tn1[l,]=Tkn1.2-Tkn1.1
}
tn1=tn1^2
tn1star=colMeans(tn1)
k1star=which.min(tn1star)

rho=secorder[2]
Exp=2/(1-2*rho)
Exp2=(-2*rho)/(1-2*rho)
k0star=floor((1-2^rho)^(Exp)*k1star*(n/n1)^(Exp2))+1
u=sort(data,decreasing=TRUE)[k0star]
ti=1/helphill(k0star,1)
list=list(sec.order.par=secorder,k0=k0star,threshold=u,tail.index=ti)
list
}
