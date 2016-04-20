fact2vec <-
function(fact,form){
#fact<- "organ"
fa <- get(fact)
lev <- levels(fa)
vec <- sapply(lev[-1],function(le) paste(fact,le,sep=""))
return(vec)
}
