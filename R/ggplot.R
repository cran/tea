#' Gerstengarbe Plot
#'
#' Performs a sequential Mann-Kendall Plot also known as Gerstengarbe Plot.
#' @param data vector of sample data
#' @param nexceed number of exceedances. Default is the minimum of the data to make sure the whole dataset is considered.
#' @details The Gerstengarbe Plot, referring to Gerstengarbe and Werner (1989), is a sequential version of the Mann-Kendall-Test. This test searches for change points within a time series. This method is adopted for finding a threshold in a POT-model. The basic idea is that the differences of order statistics of a given dataset behave different between the body and the tail of a heavy-tailed distribution. So there should be a change point if the POT-model holds. 
#' To identify this change point the sequential test is done twice, for the differences from start to the end of the dataset and vice versa. The intersection point of these two series can then be associated with the change point of the sample data. For more informations see references.
#' @return
#' \item{k0}{optimal number of upper order statistics, i.e. the change point of the dataset}
#' \item{threshold}{the corresponding threshold}
#' \item{tail.index}{the corresponding tail index}
#' @references Gerstengarbe, F.W. and Werner, P.C. (1989). A method for statistical definition of extreme-value regions and their application to meteorological time series. \emph{Zeitschrift fuer Meteorologie}, \bold{39(4)}, 224--226.
#' @references Cebrian, A., and Denuit, M. and Lambert, P. (2003). Generalized pareto fit to the society of actuaries large claims database. \emph{North American Actuarial Journal}, \bold{7(3)}, 18--36.
#' @examples
#' data(danish)
#' ggplot(danish)
#' @section Authors:
#' Ana Cebrian
#' Johannes Ossberger
#' @section Acknowledgements:
#' Great thanks to A. Cebrian for providing a basic version of this code.
#' @export

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
