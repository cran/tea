eye <-
function(data,ws=0.01,epsilon=0.3,h=0.9){

n=length(data)
w=floor(ws*n)

i=1:(n-1)
x=sort(data,decreasing=TRUE)
gamma=(cumsum(log(x[i]))/i)-log(x[i+1]) 
alpha=1/gamma

count=0;erg=c()
for (k in 2:(length(alpha)-w)){
  for (i in 1:length(w)){
    if (alpha[k+i]<(alpha[k]+epsilon) && alpha[k+i]>(alpha[k]-epsilon)) {
      count=count+1
  } else {
    count=count
  }
  }
  erg[k]=count/w
}
erg=erg>h
k0=min(which(erg==1))
u=x[k0]
ti=alpha[k0]
list=list(k0=k0,threshold=u,tail.index=ti)
list
}
