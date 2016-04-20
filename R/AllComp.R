AllComp <- function(vars){
if(is.factor(get(vars))){
fa <- get(vars)
lev <- levels(fa)
vec <- sapply(lev[-1],function(le) paste(vars,le,sep=""))
} else vec <- vars
nlev <- length(vec)
lincombvec <- c()
for(i in 1:(nlev-1)){
    for(j in (i+1):nlev){
    lc <- inla.make.lincomb(tem=-1, tem2=1)
    names(lc[[1]][[1]]) = vec[i]
    names(lc[[1]][[2]]) = vec[j]
    names(lc) <- paste(vec[j],"min",vec[i],sep="")
    lincombvec <- c(lincombvec,lc)
    }
}
return(lincombvec)
}
