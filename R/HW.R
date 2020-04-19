#' Minimizing the AMSE of the Hill estimator with respect to k
#'
#' An Implementation of the procedure proposed in Hall & Welsh (1985) for obtaining the optimal number of upper order statistics \code{k} for the Hill estimator by minimizing the AMSE-criterion.
#' @param data vector of sample data
#' @details The optimal number of upper order statistics is equivalent to the number of extreme values or, if you wish, the number of exceedances in the context of a POT-model like the generalized Pareto distribution. This number is identified by minimizing the AMSE criterion with respect to \code{k}. The optimal number, denoted \code{k0} here, can then be associated with the unknown threshold \code{u} of the GPD by choosing \code{u} as the \code{n-k0}th upper order statistic. For more information see references.
#' @return
#' \item{second.order.par}{gives an estimation of the second order parameter \code{rho}.}
#' \item{k0}{optimal number of upper order statistics, i.e. number of exceedances or data in the tail}
#' \item{threshold}{the corresponding threshold}
#' \item{tail.index}{the corresponding tail index}
#' @references Hall, P. and Welsh, A.H. (1985). Adaptive estimates of parameters of regular variation. \emph{The Annals of Statistics}, \bold{13(1)}, 331--341.
#' @examples
#' data(danish)
#' HW(danish)
#' @export
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
