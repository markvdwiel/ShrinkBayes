myinla.expectation <-
function (fun, marginal, ...) 
{  

    if(is.null(marginal)) return(NA) else {
    
    if(nrow(marginal)>200) xx <- list(x=marginal[,1],y=marginal[,2]) else 
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
    fun = match.fun(fun)
    ff = fun(xx$x[1:n]) * xx$y[1:n]
    #fff = fun(xx$x)
    #ff = fun(xx$x[1:n], ...) * xx$y[1:n]
    e = sum(sum(ff[i.0]) + 4 * sum(ff[i.4]) + 2 * sum(ff[i.2]))
    ff = 1 * xx$y[1:n]
    e.1 = sum(sum(ff[i.0]) + 4 * sum(ff[i.4]) + 2 * sum(ff[i.2]))
    return(e/e.1)
    }
}
