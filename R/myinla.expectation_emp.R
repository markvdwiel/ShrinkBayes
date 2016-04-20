myinla.expectation_emp <-
function(dens, marginal,lengthmin=200) 
{
    if(nrow(marginal)>=lengthmin) xx <- list(x=marginal[,1],y=marginal[,2]) else 
    {
        xx = try(myinla.smarginal(marginal),silent=TRUE)
        if(class(xx)=="try-error") xx <- list(x=marginal[,1],y=marginal[,2])
        if(class(xx)!="try-error")  if(max(xx$y)>10) xx <- list(x=marginal[,1],y=marginal[,2])
    }
    n = length(xx$x)
    if (n%%2 == 0) 
        n = n - 1
    i.0 = c(1, n)
    i.4 = seq(2, n - 1, by = 2)
    i.2 = seq(3, n - 2, by = 2)
    
    #ff = fun(xx$x[1:n], ...) * xx$y[1:n]
    nx <- dens[,1]
    ny <- dens[,2]
    minx <- nx[1]
    maxx <- nx[length(nx)]
    ndx <- function(ex) {if(ex <= minx | ex >= maxx) {return(0)} else
        {
        #ex<-0
        wh <- min(which(nx>=ex))-1
        res <- ny[wh] + ((ny[wh+1]-ny[wh])/(nx[wh+1]-nx[wh]+10^(-10)))*(ex-nx[wh])
        return(res)
        }
    }
    ff = sapply(xx$x[1:n],ndx) * xx$y[1:n]
    e = sum(sum(ff[i.0]) + 4 * sum(ff[i.4]) + 2 * sum(ff[i.2]))
    ff = 1 * xx$y[1:n]
    e.1 = sum(sum(ff[i.0]) + 4 * sum(ff[i.4]) + 2 * sum(ff[i.2]))
    return(e/e.1)
}
