myinla.qmarginal <-
function (p, marginal, len = 1024) 
{ 
    if(nrow(marginal)>200) {
    lin <- TRUE
    marginal <- list(x=marginal[,1],y=marginal[,2])
    f = linear.spline(marginal)
    } else { #apply spline
    marginal <- marginal[marginal[,2]>=10^{-20},]
    lin <- FALSE
    f = INLA:::inla.sfmarginal(myinla.smarginal(marginal))
    }
    xx = seq(f$range[1], f$range[2], length = len)
    d = cumsum(exp(f$fun(xx)))
    d = d/d[length(d)]
    
    if(lin) {fq = approxfun(d,xx)} else {fq = splinefun(d,xx, method = "monoH.FC")}
   
    n = length(p)
    pp = pmin(pmax(p, rep(0, n)), rep(1, n))
    
    toret <- try(fq(pp),silent=TRUE)
    if(class(toret)=="try-error") toret <- sapply(pp,function(ppi) xx[which.max(d[d<=ppi])])
    return(toret)
}
