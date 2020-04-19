#' Sample Path Stability Algorithm
#'
#' An Implementation of the heuristic algorithm for choosing the optimal sample fraction proposed in Caeiro & Gomes (2016), among others.
#' @param data vector of sample data
#' @param j digits to round to. Should be \code{0} or \code{1} (default)
#' @details The algorithm searches for a stable region of the sample path, i.e. the plot of a tail index estimator with respect to \code{k}. This is done in two steps. First the estimation of the tail index for every \code{k} is rounded to \code{j} digits and the longest set of equal consecutive values is chosen. For this set the estimates are rounded to \cite{j+2} digits and the mode of this subset is determined. The corresponding biggest k-value, denoted \code{k0} here, is the optimal number of data in the tail. 
#' @return
#' \item{k0}{optimal number of upper order statistics, i.e. number of exceedances or data in the tail}
#' \item{threshold}{the corresponding threshold}
#' \item{tail.index}{the corresponding tail index}
#' @references Caeiro, J. and Gomes, M.I. (2016). Threshold selection in extreme value analysis. \emph{Extreme Value Modeling and Risk Analysis:Methids and Applications}, 69--86.
#' @references Gomes, M.I. and Henriques-Rodrigues, L. and Fraga Alves, M.I. and Manjunath, B. (2013). Adaptive PORT-MVRB estimation: an empirical comparison of two heuristic algorithms. \emph{Journal of Statistical Computation and Simulation}, \bold{83}, 1129--1144.
#' @references Gomes, M.I. and Henriques-Rodrigues, L. and Miranda, M.C. (2011). Reduced-bias location-invariant extreme value index estimation: a simulation study. \emph{Communications in Statistic-Simulation and Computation}, \bold{40}, 424--447.
#' @examples
#' data(danish)
#' PS(danish)
#' @export
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
