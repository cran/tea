PS <-
function(data,j=1){
n=length(data)
x1=sort(data, decreasing=TRUE)
i=1:(n-1)
h=(cumsum(log(x1[i]))/i)-log(x1[i+1])

help=round(h,j)
l=rle(help)$lengths
Max=which(l==max(l))
kmin=sum(l[1:(Max-1)])+1
kmax=kmin+max(l)-1
run=kmax-kmin

K=round(h[kmin:kmax],j+2)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
if (sum(duplicated(K))==0) {
  modus=K[length(K)]
} else{
  modus=Mode(K)
}

k0star=max(which(round(h,j+2)==modus))
u=x1[k0star]
ti=1/h[k0star]
list=list(k0=k0star,threshold=u,tail.index=ti)
list
}
