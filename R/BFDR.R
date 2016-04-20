BFDR <-  function(lfdr1,lfdr2=NULL,multcomp=FALSE){
if(multcomp){
    if(class(lfdr1)=="numeric") {lfdr1 <- matrix(lfdr1,ncol=1);if(!is.null(lfdr2)){lfdr2 <- matrix(lfdr2,ncol=1)}}
    if(!is.null(lfdr2)) lfdr <- apply(cbind(lfdr2,lfdr1),1,min) else lfdr <- apply(lfdr1,1,min)
    whna <- which(is.na(lfdr))
    if(length(whna)>0) lfdr[whna]<-1
    allthr <- sort(lfdr,index.return=T)
    BFDR0 <- rep(0,length(allthr$x))
    BFDR0[allthr$ix] <-cumsum(allthr$x)/(1:length(allthr$x))
    BFDR<- BFDR0
    return(BFDR)
} else {
    if(is.null(lfdr2)){ #one-sided testing
        if(class(lfdr1)=="numeric") lfdr1 <- matrix(lfdr1,ncol=1)
        best <- lfdr1
        whna <- which(is.na(best))
        if(length(whna)>0) best[whna]<-1
        allthr <- sort(best,index.return=T)
        BFDR <- rep(0,length(allthr$x))
        BFDR[allthr$ix] <-cumsum(allthr$x)/(1:length(allthr$x))
        BFDRres0 <- matrix(BFDR,ncol=ncol(lfdr1))
        colnames(BFDRres0) <- colnames(lfdr1)
        BFDRres <- matrix(nrow=length(best),ncol=ncol(BFDRres0))
        BFDRres<-BFDRres0
        return(BFDRres)
    } else { #two-sided
        if(class(lfdr1)=="numeric") {
        best <- mapply(min,lfdr2,lfdr1)
        lfdr1 <- matrix(lfdr1,ncol=1);lfdr2 <- matrix(lfdr2,ncol=1)
        } else {
        best <- mapply(min,lfdr2,lfdr1)
        }
        whna <- which(is.na(best))
        if(length(whna)>0) best[whna]<-1
        allthr <- sort(best,index.return=T)
        BFDR <- rep(0,length(allthr$x))
        BFDR[allthr$ix] <-cumsum(allthr$x)/(1:length(allthr$x))
        BFDRres0 <- matrix(BFDR,ncol=ncol(lfdr1))
        colnames(BFDRres0) <- colnames(lfdr1)
        BFDRres <- matrix(nrow=length(best),ncol=ncol(BFDRres0))
        BFDRres<-BFDRres0
        return(BFDRres)
    }
}
}
