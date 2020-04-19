#' Adaptive choice of the optimal sample fraction in tail index estimation
#'
#' An implementation of the minimization criterion proposed in Reiss & Thomas (2007).
#' @param data vector of sample data
#' @param beta a factor for weighting the expression below. Default is set to \code{beta=0}
#' @param kmin gives a minimum value for \code{k}. Default ist set to \code{kmin=2}
#' @details The procedure proposed in Reiss & Thomas (2007) chooses the lowest upper order statistic \code{k} to minimize the expression
#' \code{1/k sum_i=1^k i^beta |gamma_i-median(gamma_1,...,gamma_k)|}
#' or an alternative of that by replacing the absolute deviation with a squared deviation and the median just with \code{gamma_k}, where \code{gamma} denotes the Hill estimator
#' @return
#' \item{k0}{optimal number of upper order statistics, i.e. number of exceedances or data in the tail for both metrics, i.e. the absolute and squared deviation.}
#' \item{threshold}{the corresponding thresholds.}
#' \item{tail.index}{the corresponding tail indices}
#' @references Reiss, R.-D. and Thomas, M. (2007). Statistical Analysis of Extreme Values: With Applications to Insurance, Finance, Hydrology and Other Fields. \emph{Birkhauser, Boston}.
#' @examples
#' data(danish)
#' RT(danish)
#' @export
RT <-
function(data,beta=0,kmin=2){
n=length(data)

x1=sort(data, decreasing=TRUE)
i=1:(n-1)
j=i^beta
h=(cumsum(log(x1[i]))/i)-log(x1[i+1])

med=c();
for (i in 1:(n-1)){
  med[i]=median(h[1:i])
}
erg1=c(); erg2=c()
for (k in 1:(n-1)){
erg1[k]=sum(j[1:k]*abs(h[1:k]-med[k]))/k
erg2[k]=sum(j[1:k]*(h[1:k]-h[k])^2)/k
}
rt1=which.min(erg1[kmin:(n-1)])
rt2=which.min(erg2[kmin:(n-1)])
list=list(k0=c(rt1,rt2),threshold=c(x1[rt1],x1[rt2]),tail.index=c(1/h[rt1],1/h[rt2]))
list
}
