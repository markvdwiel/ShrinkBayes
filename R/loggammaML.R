loggammaML <-
function(xx, startval=c(1,1)){
    xx<-xx[xx!=Inf]
    loglikloggam <- function(a=1,b=1){sum(-dlgamma(xx, location=0, 1/b, a, log=T))}  #a is shape, scale = 1/b
    res1 <- mle(loglikloggam,start=list(a=startval[1],b=startval[2]),lower=c(0.1,0.1),method="L-BFGS-B")
    if(mode(res1)=="S4") coefs <- res1@coef else coefs <- coef(res1)
    return(c(coefs[1],coefs[2])) #shape and rate
}
