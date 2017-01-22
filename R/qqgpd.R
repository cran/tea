qqgpd<-
function(data,nextremes,scale,shape){ 
  xstat=sort(data,decreasing=T)[1:nextremes]
  u=sort(data,decreasing=T)[nextremes]
  i=ppoints(nextremes)
  x=qgpd(i,loc=u,scale=scale,shape=shape) 
  
  sim.new.quant=function(l){
    simdata=rgpd(nextremes,loc=u,scale=scale,shape=shape)
    sort(simdata)
  }
  sim.q=sapply(1:1000,sim.new.quant)
  ci.q=apply(sim.q,quantile,MARGIN=1,probs=c(0.025,0.975))
    
  axislims=c(u,max(xstat))
  plot(x,sort(xstat),pch="x",xlab="theoretical quantiles",ylab="empirical quantiles",
       xlim=axislims,ylim=axislims)
  abline(0,1,col="red")
  lines(x,ci.q[1,],lty=2)
  lines(x,ci.q[2,],lty=2)
}