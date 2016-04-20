myinla.pmarginal <-
function(x,marginal,log=FALSE, len = 1024){
    if(nrow(marginal)<=200) {
        postspl <- try(myinla.smarginal(marginal))
        if(class(postspl)=="try-error") postspl <- list(x=marginal[,1],y=marginal[,2])
        if(class(postspl)!="try-error")  if(max(postspl$y)>10) postspl <- list(x=marginal[,1],y=marginal[,2])
        f = INLA:::inla.sfmarginal(myinla.smarginal(marginal))
        xx = seq(f$range[1], f$range[2], length = nrow(marginal))
        d = cumsum(exp(f$fun(xx)))
        d = d/d[length(d)]
        fq = splinefun(xx, d, method = "monoH.FC")
        n = length(x)
        xx = pmin(pmax(x, rep(f$range[1], n)), rep(f$range[2], n))
        return(fq(xx)) 
    } else { #use linear splines when splines have already been used    
        postspl <- list(x=marginal[,1],y=marginal[,2]) 
        f = linear.spline(postspl)
        xx = seq(f$range[1], f$range[2], length = nrow(marginal))
        d = cumsum(exp(f$fun(xx)))
        d = d/d[length(d)]
        fq = approxfun(xx, d)
        n = length(x)
        xx = pmin(pmax(x, rep(f$range[1], n)), rep(f$range[2], n))
        return(fq(xx))  
    }  
}
