#' Minimizing the AMSE of the Hill estimator with respect to k
#'
#' Gives the optimal number of upper order statistics \code{k} for the Hill estimator by minimizing the AMSE-criterion.
#' @param data vector of sample data
#' @details The optimal number of upper order statistics is equivalent to the number of extreme values or, if you wish, the number of exceedances in the context of a POT-model like the generalized Pareto distribution. This number is identified by minimizing the AMSE criterion with respect to \code{k}. The optimal number, denoted \code{k0} here, can then be associated with the unknown threshold \code{u} of the GPD by choosing \code{u} as the \code{n-k0}th upper order statistic. For more information see references.
#' @return 
#' \item{second.order.par}{gives an estimation of the second order parameter \code{beta} and \code{rho}.}
#' \item{k0}{optimal number of upper order statistics, i.e. number of exceedances or data in the tail}
#' \item{threshold}{the corresponding threshold}
#' \item{tail.index}{the corresponding tail index}
#' @references Caeiro, J. and Gomes, M.I. (2016). Threshold selection in extreme value analysis. \emph{Extreme Value Modeling and Risk Analysis:Methids and Applications}, 69--86.
#' @examples
#' data(danish)
#' dAMSE(danish)
#' @export
dAMSE <-
function(data){
n=length(data)
k=c(floor(n^0.995),floor(n^0.999))

helphill=function (k,j=1) {
  xstat = sort(data, decreasing = TRUE)
  xihat = mean((log(xstat[1:k]) - log(xstat[k + 1]))^j)
  xihat
}
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

Exp=1/(1-2*rho)
Z=(1-rho)^2*n^(-2*rho)
N=-2*rho*beta^2
k0=floor((Z/N)^Exp)
u=xstat[k0]
ti=1/helphill(k0)
list=list(second.order.par=c(beta,rho),k0=k0,threshold=u,tail.index=ti)

list
}
