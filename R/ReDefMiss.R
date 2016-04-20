ReDefMiss <- function(para,NAzero=F){
para <- get(para)
if(is.factor(para)) {
parach <- as.character(para) 
parach[is.na(parach)] <- "0"
para <- factor(parach)
} else 
    {
    if(NAzero) para[is.na(para)] <- 0 else {
        if(sum(is.na(para)>0)){
        if(length(unique(para))==3) {  #binary coding,NA counts as one level
        minval <- min(para,na.rm=T);maxval <- max(para,na.rm=T)
        para[para==minval] <- -1
        para[para==maxval] <- 1
        para[is.na(para)] <- 0
        } else {  #continuous coding, center first
        para <- para - median(para,na.rm=T)
        para[is.na(para)] <- 0
        }
        }
    }
}
return(para)
}
