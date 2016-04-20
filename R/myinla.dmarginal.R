myinla.dmarginal <-
function(x,marginal,log=FALSE,thrnrow=200){
#x<-marg1[,1];marginal<-marg2;thrnrow=140
    if(nrow(marginal)>thrnrow) {
        postspl <- list(x=marginal[,1],y=marginal[,2]) 
        f = try(linear.spline(postspl))
        } else {
        postspl <- try(myinla.smarginal(marginal))
        if(class(postspl)=="try-error") postspl <- list(x=marginal[,1],y=marginal[,2])
        if(class(postspl)!="try-error")  if(max(postspl$y)>10) postspl <- list(x=marginal[,1],y=marginal[,2])
        f = try(INLA:::inla.sfmarginal(postspl))
        }
    if(class(f)=="try-error") return(rep(0,length(x))) else {
    n = length(x)
    d = numeric(n)
    for (i in 1:n) {
        if (x[i] >= f$range[1] && x[i] <= f$range[2]) {
            if (log) 
                d[i] = f$fun(x[i])
            else d[i] = exp(f$fun(x[i]))
        }
        else {
            if (log) 
                d[i] = -Inf
            else d[i] = 0
        }
    }
    return(d)
    }
}
