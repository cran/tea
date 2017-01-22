TH<-
function(data,thresholds){
u=thresholds

pearson.test=function (x, n.classes = ceiling(2 * (n^(2/5))), adjust = TRUE) 
{
  DNAME <- deparse(substitute(x))
  x <- x[complete.cases(x)]
  n <- length(x)
  if (adjust) {
    dfd <- 2
  }
  else {
    dfd <- 0
  }
  num <- floor(1 + n.classes * pnorm(x, mean(x), sd(x)))
  count <- tabulate(num, n.classes)
  prob <- rep(1/n.classes, n.classes)
  xpec <- n * prob
  h <- ((count - xpec)^2)/xpec
  P <- sum(h)
  pvalue <- pchisq(P, n.classes - dfd - 1, lower.tail = FALSE)
  RVAL <- list(statistic = c(P = P), p.value = pvalue, method = "Pearson chi-square normality test", 
               data.name = DNAME, n.classes = n.classes, df = n.classes - 
                 1 - dfd)
  class(RVAL) <- "htest"
  return(RVAL)
}

shape=c();scale=c();exc=c()
for (l in 1:length(u)) {
  est=gpdFit(data,threshold=u[l])
  exc[l]=est$n.exceed
  shape[l]=est[[4]][2,1]
  scale[l]=est[[4]][1,1]
  rm(est)
}

tu=c();
for (j in 1:length(u)) {
  tu[j]=scale[j]-shape[j]*u[j] 
}
t=diff(tu)     

p.v=c();
for (q in 1:(length(t)-2)){
  p.v[q]=pearson.test(t[q:length(t)])$p.value
}
stops=pSeqStop(p.v)
names=c("testnum","threshold","num.above","p.values",
        "ForwardStop","StrongStop","est.scale","est.shape")
erg=matrix(ncol=8,nrow=length(u)-3)
colnames(erg)=names
erg[,1]=1:(length(u)-3)
erg[,2]=u[1:(length(u)-3)]
erg[,3]=exc[1:(length(u)-3)]
erg[,4]=c(p.v)
erg[,5]=c(stops$ForwardStop)
erg[,6]=c(stops$StrongStop)
erg[,7]=scale[1:(length(u)-3)]
erg[,8]=shape[1:(length(u)-3)]

erg
}