mindist<-
function(data,ts=0.15,method="mad"){
  xstat=sort(data,decreasing=TRUE)
  n=length(data)
  T=floor(n*ts)
  i=1:(n-1)
  h=(cumsum(log(xstat[i]))/i)-log(xstat[i+1]) 
  xstat=sort(data)
  A=matrix(ncol=T-1,nrow=T-1)
  for (k in 1:(T-1)){
    for (j in 1:(T-1)){
      A[k,j]=abs( (((k/j)*xstat[n-k+1]^(1/h[k]))^h[k]) - xstat[n-j] )
    }
  }
  if (method=="mad"){
    M=rowMeans(A)
    kstar=which.min(M)
    u=rev(xstat)[kstar]
    list=list(k0=kstar,threshold=u,tail.index=1/h[kstar])
    return(list)   
  } 
  if (method=="ks"){
      rowMax <- function (rowData) {
      apply(rowData, MARGIN=c(1), max)
    }
    M=rowMax(A)
    kstar=which.min(M)
    u=rev(xstat)[kstar]
    list=list(k0=kstar,threshold=u,tail.index=1/h[kstar])
    return(list)
  }
  if (method!="mad" && method!="ks"){
    warning("method should be one of 'mad' or 'ks'")
  }  
}