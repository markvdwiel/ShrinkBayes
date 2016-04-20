BaselineDef <-
function(fact,baselinegroup=as.character(levels(get(fact)))[1])
{
fa <- as.vector(get(fact))
wh <- which(fa==baselinegroup)
if(length(wh)==0){
print("Baseline level not found, factor levels remain unchanged")
} else {
wh0 <- which(fa=="0")
if(length(wh0)>0){
print("Level 0 is already defined. This is now redefined to level 0b")
fa[wh0] <- "0b"
}
fa[wh] <- 0
}
return(as.factor(fa))
}
