gammaML <-
function(xx, startval=c(1,1)){
    xx<-xx[xx>0]
    loglikgam <- function(a=1,b=1){sum(-dgamma(xx, shape=a, scale=1/b, log=T))}  #k is shape, scale = 1/rate
    res1 <- mle(loglikgam,start=list(a=startval[1],b=startval[2]),lower=c(0.1,0.1),method="L-BFGS-B")
    if(mode(res1)=="S4") coefs <- res1@coef else coefs <- coef(res1)
    return(c(coefs[1],coefs[2])) 
}

