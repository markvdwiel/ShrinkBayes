mixture2loggammaML <-
function(x,thrfit = 8){
    p0 <- sum(abs(x)==Inf)/length(x)
    x <- x[x!=Inf]
    loglikloggam <- function(a=1,b=1){sum(-dlgamma(xx, location=0, 1/b, a, log=T))}  #k is shape, scale = 1/rate
    res1 <- mle(loglikloggam,start=list(a=1,b=1),lower=c(0.1,0.1),method="L-BFGS-B")
    if(mode(res1)=="S4") coefs <- res1@coef else coefs <- coef(res1)
    return(c(p0,coefs[1],coefs[2])) #p0,shape and rate
}
