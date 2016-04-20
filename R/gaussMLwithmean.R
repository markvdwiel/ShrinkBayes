gaussMLwithmean <-
function(x,startval=NULL,thrfit=5){ 
    x <- x[abs(x)<=thrfit]
    out <-c(mn=mean(x),stdv=1/(sd(x))^2) #meanprec
    return(out)
}

