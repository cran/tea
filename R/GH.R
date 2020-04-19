#' A Bias-based procedure for Choosing the Optimal Threshold
#'
#' An Implementation of the procedure proposed in Guillou & Hall(2001) for selecting the optimal threshold in extreme value analysis.
#' @param data vector of sample data
#' @details The procedure proposed in Guillou & Hall (2001) is based on bias reduction. Due to the fact that the log-spacings of the order statistics are approximately exponentially distributed if the tail of the underlying distribution follows a Pareto distribution, an auxilliary statistic with respect to \code{k} is implemented with the same properties. The method then behaves like an asymptotic test for mean \code{0}. If some critical value \code{crit} is exceeded the hypothesis of zero mean is rejected. Thus the bias has become too large and the assumed exponentiality and therefore the assumed Pareto tail can not be hold.  From this an optimal number of \code{k} can be found such that the critical value is not exceeded. This optimal number, denoted \code{k0} here, is equivalent to the number of extreme values or, if you wish, the number of exceedances in the context of a POT-model like the generalized Pareto distribution. \code{k0} can then be associated with the unknown threshold \code{u} of the GPD by 
#' coosing \code{u} as the \code{n-k0}th upper order statistic. For more information see references.
#' @return
#' \item{k0}{optimal number of upper order statistics, i.e. number of exceedances or data in the tail}
#' \item{threshold}{the corresponding threshold}
#' \item{tail.index}{the corresponding tail index}
#' @references Guillou, A. and Hall, P. (2001). A Diagnostic for Selecting the Threshold in Extreme Value Analysis. \emph{Journal of the Royal Statistical Society}, \bold{63(2)}, 293--305.
#' @examples
#' data(danish)
#' GH(danish)
#' @export
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
