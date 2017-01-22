ggplot <-
function(data,nexceed=min(data)-1)
{
	      nmin <- NULL
        n<-NULL
        np<-NULL
        u<-NULL
        up<-NULL
        yp<-NULL
        exc<-data[data>nexceed]
#exc contien los datos de la serie original superiores al umbral de referencia
#como vamos a trabajar con diferencias de estos datos, las diferencias son invariantes
#respecto al umbral y no tenemos que restar el umbral  a cada observacion original
#el valor de l adiferncia sera el mismo
        sexc<--sort(-exc)
#sexc contien los datos de la serie original superiores al umbral de referencia
#ordenados de mayor a menor

        y<--diff(sexc)
        ll<-length(y)
#ll es la longitud de y



	yp<-y[ll:1]
#yp contiene la serie de diferencias en sentido inverso
        j<-1:ll
	n<-NULL
	npi<-NULL
        for(k in 1:ll) 
	{
                n[k]<-sum(as.numeric((j<k)&(y<y[k])))
                npi[k]<-sum(as.numeric((j<k)&(yp<yp[k])))

 	}
#n  sucht den Rang von y[k] innerhalb der ersten k Differenzen
#npi  analog
        t<-cumsum(n)
        tpi<-cumsum(npi)	
#t contiene en cada fila i el estadistico ti, suma de los n(k) anteriores
#tpi contiene en cada fila i el estadistico tpi, suma de los npi(k) anteriores
# es decir el estadistico ti calculado con la serie de diferencias invertida

        i<-1:ll
        u<-((t-(i*(i-1))/4)*(72**0.5))/((i*(i-1)*(2*i+5))**0.5)
        upi<--((tpi-(i*(i-1))/4)*(72**0.5))/((i*(i-1)*(2*i+5))**0.5)
#en u estan los estadisticos ti estandarizados
#en upi estan los estadisticos tpi estandarizados
# por la definicion de estos hay que invertir el vector, que se almacena en up
# u y up se representan frente a i
	u[1]<-0
	upi[1]<-0

	up<-upi[length(upi):1]

        maxig<-max(u,up)
	minig<-min(u,up)

	      plot(i, u, type = "l", ylim=c(minig,maxig),xlab = "k", ylab = "Up, Ur")
	#,main=paste('Gerstengarbe plot')) 
	par(lty=2)
        lines(i,up,type="l")
	par(lty=1)

#Schnittpunkte berechnen

	difuup<-as.numeric((u-up)>0)
	pc<-i[c(0,diff(difuup))!=0]
	
	fu<-function()
	{
	plot(i, u, type = "l", ylim=c(-7,7),xlab = "(i)", ylab = "u, up")
	par(lty=2)
        lines(i,up,type="l")
	par(lty=1)	
	}

	pv<-2*(1-pnorm(abs(u[pc]),0,1))

  helphill=function (k) {
    xstat = sort(data, decreasing = TRUE)
    xihat = mean((log(xstat[1:k]) - log(xstat[k + 1])))
    xihat
  }
 
  ti=1/sapply(pc,helphill)

  list=list(k0=pc,p.values=pv,threshold=sexc[pc],tail.index=ti)
  list
		}
