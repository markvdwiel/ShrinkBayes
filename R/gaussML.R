gaussML <-
function(x,startval=NULL,thrfit=5){ 
    x <- x[abs(x)<=thrfit]
    loglikgauss <- function(stdv=1){sum(-dnorm(x,mean=0, sd=stdv, log=T))} 
    if(is.null(startval)) res1 <- mle(loglikgauss,lower=c(0.01),method="L-BFGS-B") else res1 <- mle(loglikgauss,start=list(stdev=startval),lower=c(0.01),method="L-BFGS-B")
    if(mode(res1)=="S4") coefs <- res1@coef else coefs <- coef(res1)
    out <- c(0,1/coefs[1]^2) #meanprec
    return(out)
}

