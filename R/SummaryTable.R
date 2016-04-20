SummaryTable <- function(poster,BFDRthr=0.1,diffthr = 0,direction="two-sided",pointmass=0,ndigit=3,ncpus=1){  #direction=c("two-sided","greater", "smaller", "equal")
#poster <- posteriors;diffthr = 0; direction = "equal";ncpus=1;updateby=10000;pointmass=0;BFDRthr=1
lfdrs <- SummaryWrap(poster,thr=diffthr,direction=direction,ncpus=ncpus)
BFDRs <- BFDR(lfdrs)
postmean <- SummaryWrap(poster, summary="postmean",ncpus=ncpus)
postpi0 <- SummaryWrap(poster,summary="postpi0",ncpus=ncpus)
namsfdr <- colnames(lfdrs)
nams <- colnames(postmean)
colnames(lfdrs) <- sapply(namsfdr,function(nm) paste("lfdr","_",nm,sep=""))
colnames(BFDRs) <- sapply(namsfdr,function(nm) paste("BFDR","_",nm,sep=""))
if(!is.null(postmean)) colnames(postmean) <- sapply(nams,function(nm) paste("postmean","_",nm,sep=""))
if(ncol(postpi0)>1) colnames(postpi0) <- sapply(nams,function(nm) paste("prob0","_",nm,sep=""))
if(ncol(BFDRs)>1) {
BFDRmin <- apply(BFDRs,1,min)
wh <- which(BFDRmin <= BFDRthr)
} else wh <- which(BFDRs <= BFDRthr)
if(!is.null(postmean)) ST <- data.frame(lfdrs,BFDRs,postpi0,postmean)[wh,] else ST <- data.frame(lfdrs,BFDRs,postpi0)[wh,]
ST <- round(ST,3)
ST <- data.frame(index=wh,ST)
return(ST)
}
