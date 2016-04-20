SummaryWrap <-
function(posteriorlist, thr=0,ncpus=1,updateby=10000,summary="lfdr",direction = "two-sided",pointmass=0){
#direction = ifelse(thr <= 0, "greater",   "lesser")
#posteriorlist=postshr; thr = 0; direction = "equal";ncpus=1;updateby=10000;summary="lfdr";pointmass=0
#posteriorlist <- mixtpostshr;summary="postpi0";thr = 0; direction = "equal";ncpus=1;updateby=10000;pointmass=0
#posteriorlist <- nppostshr;summary="postmean";thr = 0; direction = "equal";ncpus=1;updateby=10000;pointmass=0
posteriorlist100 <- posteriorlist[1:min(100,length(posteriorlist))]
whnotnull <- which(!is.null(posteriorlist))
ncomp <- length(posteriorlist100[[whnotnull[1]]][[1]]) 
nams <- names(posteriorlist100[[whnotnull[1]]]$postbetanon0)
absthr <- abs(thr)

if((is.null(posteriorlist[[1]]$postbeta0) | posteriorlist[[1]]$postbeta0[1] == 0) & thr == 0) {
print("Advise is to use a threshold thr > 0 when a prior without a spike on 0 is used")
}

if(summary=="postmean" & is.null(posteriorlist[[1]]$postbetanon0)){
print("Not possible to compute posterior means, returning NULL")
return(NULL)
}

if(is.null(posteriorlist[[1]]$postbeta0)){
    summarylist <- function(poster,thr,ncomp){
    if(is.null(poster)) return(rep(NA,ncomp)) else {
    poster1 <- poster[[1]]
    if(summary=="lfdr") {
    if(direction=="greater" | direction=="lesser"){
    lfdr1 <- sapply(poster1,function(x) lfdr(x,thr=thr,direction=direction))   
    }
    if(direction=="two-sided"){
    lfdr1a <- (1-poster0)*sapply(poster1,function(x) lfdr(x,thr=absthr,direction="lesser"))
    lfdr1b <- (1-poster0)*sapply(poster1,function(x) lfdr(x,thr=-absthr,direction="greater"))
    lfdr1 <- mapply(min, lfdr1a, lfdr1b)
    }
    return(lfdr1)
    }
    if(summary=="postmean") return(sapply(poster1,function(x) postmean(x)))
    if(summary=="postpi0") return(sapply(poster1,function(x) 0))
    }
    } 
} else {
summarylist <- function(poster,thr,ncomp){
    if(is.null(poster)) return(rep(NA,ncomp)) else {
    poster1 <- poster$postbetanon0
    poster0 <- poster$postbeta0
    if(summary=="lfdr") {
        if(direction=="equal" | poster0==1) return(poster0) else {
            if(direction=="greater" | direction=="lesser"){
            lfdr1 <- (1-poster0)*sapply(poster1,function(x) lfdr(x,thr=thr,direction=direction))
            if((thr<=pointmass & direction =="greater") | (thr>=pointmass & direction =="lesser")) lfdr1 <- lfdr1 + poster0
            }
            if(direction=="two-sided"){
            lfdr1a <- (1-poster0)*sapply(poster1,function(x) lfdr(x,thr=absthr,direction="lesser")) + poster0
            lfdr1b <- (1-poster0)*sapply(poster1,function(x) lfdr(x,thr=-absthr,direction="greater")) + poster0
            lfdr1 <- mapply(min, lfdr1a, lfdr1b)
            }
            return(lfdr1)
            }
            }
            
    if(summary=="postmean") {if(poster0==1) return(pointmass) else return(pointmass*poster0 + (1-poster0)*sapply(poster1,function(x) postmean(x)))}
    if(summary=="postpi0") return(poster0)  
    } 
}
}


if(ncpus==1) {if(ncomp>1) res <- t(sapply(posteriorlist,summarylist,thr=thr,ncomp=ncomp)) else 
res <- matrix(sapply(posteriorlist,summarylist,thr=thr,ncomp=ncomp),ncol=1)} else
{
sfInit(parallel=TRUE,cpus=ncpus)
sfLibrary(INLA) 
el <- length(posteriorlist)
if(el<=updateby){
    print("Exporting to slaves")
    #sfExport("posteriorlist","lfdr","thr","ncomp","summarylist","myinla.pmarginal","myinla.smarginal","linear.spline")
    mysfExport(forceexport="ds")
    print("Exporting done")
    if(ncomp>1) res <- t(sfSapply(posteriorlist,summarylist,thr=thr,ncomp=ncomp)) else
    res <- matrix(sfSapply(posteriorlist,summarylist,thr=thr,ncomp=ncomp),ncol=1)
    sfRemoveAll()
} else {
print(paste("Large objects...Job cut in pieces of",updateby))
        res <- matrix(NA,nrow=el,ncol=ncomp)
        sec <- seq(1,el,by=updateby)
        secl <- length(sec)
        if(sec[secl] > el-2) sec[secl] <- el+1 else sec <- c(sec,el+1)
        secl <- length(sec)
        for(k in 1:(secl-1)){
        #k<-1
        posteriorlistk <- posteriorlist[(sec[k]:(sec[k+1]-1))]
        fnhelp <- function(ds){
        print("Exporting to slaves")
        #sfExport("ds","lfdr","thr","ncomp","summarylist","myinla.pmarginal","myinla.smarginal","linear.spline")
        mysfExport(forceexport="ds")
        print("Exporting done")
        if(ncomp>1) resk <- t(sfSapply(ds,summarylist,thr=thr,ncomp=ncomp)) else
        resk <- matrix(sfSapply(ds,summarylist,thr=thr,ncomp=ncomp),ncol=1)
        return(resk)
        }
        resk <- fnhelp(posteriorlistk)
        print(paste(sec[k+1]-1,"Cases done"))
        res[(sec[k]):(sec[k+1]-1),] <- resk
        sfRemoveAll()
        }
}
sfStop()
}
colnames(res) <- nams
if(summary=="postpi0" | summary=="lfdr") {  #collapse to one column when possible
nr <- nrow(res)
tkrow <- min(500,nrow(res))
ressel <- res[sample(1:nr,tkrow),,drop=FALSE]
uniq <- TRUE
ro <- 1
while(uniq & ro <= tkrow) {if(length(unique(ressel[ro,]))>1) {uniq <- FALSE} else {ro <- ro+1}}
if(uniq) {
    res <- res[,1,drop=FALSE]
    if(summary=="postpi0") colnames(res) <- "pi0NullModel" 
    }
}
return(res)
}
