sumplot=function (data, kmin = 5) 
{
  n = length(data)
  i = 1:(n - 1)
  xstat = sort(data, decreasing = TRUE)
  h = (cumsum(log(xstat[i]))/i) - log(xstat[i + 1])
  plot(i[kmin:(n - 1)], i[kmin:(n-1)]/h[kmin:(n - 1)], type = "l", xlab = "Order Statistics", 
       ylab = "Tail Index")
}