SpikeSlabUpdatePosterior <-
function(fitall,fitall0, p0est, modus="fixed", shrinkpara=NULL, shrinklc=NULL, ncpus=1,only0=FALSE){
continue <- TRUE
if(!is.null(fitall0)){if(!fitall0[[2]]$finalprior) continue <- FALSE}
if(!fitall[[2]]$finalprior) continue <- FALSE
if(!continue){
print("ERROR: Use finalprior = TRUE in function FitAllShrink before running SpikeSlabUpdatePosterior. Computation is aborted")
return(NULL)
} else {
#fitall=fitg;fitall0=fitg0; p0est=P0ssExact; modus="fixed"; shrinkpara="group"; shrinklc=NULL; ncpus=1;
fitall <- fitall[[1]]
fitall0 <- fitall0[[1]]
mlik <- mliks(fitall)
mlik0 <- mliks(fitall0)

if(only0) modus <- "only0"

if(modus=="fixed") {if(!is.null(shrinkpara)){if(is.factor(try(get(shrinkpara),silent=TRUE))) shrinkpara <- fact2vec(shrinkpara)}}

if(ncpus > 1) {
    sfInit(parallel=TRUE,cpus=ncpus) 
    sfLibrary(VGAM)
    sfLibrary(ShrinkBayes)
    }
    
wh <- which.min(unlist(lapply(fitall,is.null)))
callmode <- fitall[[wh]]$call


if(is.null(callmode)) { #input is NOT AN inla output object
pxbeta <- fitall  
} else {
    if(modus=="fixed") { #input is inla output object
pxbeta <- lapply(fitall, function(ex2) {
#ex2 <- fitall[[1]]
                         mhs <- ex2$marginals.fixed
                        if(!is.null(mhs)){
                            nms <- names(mhs) 
                            wh <- which(nms %in% as.vector(shrinkpara))         
                            mhswh <- mhs[wh]
                            if(!is.null(shrinklc))  {
                            mhslc <- ex2$marginals.lincomb.derived
                            nmslc <- names(mhslc)
                            whlc <- which(nmslc %in% shrinklc)
                            mhswh <- c(mhswh,mhslc[whlc])
                            }
                            whnull <- sapply(mhswh,is.null)
                            if(sum(whnull)>0) mhswh <- NULL
                            } else {
                            mhswh <- NULL
                            }
                        return(mhswh)
                        })
}


    
    if(modus=="logdisp"){  
        
        pxbeta <- lapply(fitall, function(ex2) {
                                #ex2 <- fitall[[1]]
                                mhs <- ex2$internal.marginals.hyperpar$"log size"
                                return(mhs)
                                })
    }
    
    if(modus=="random"){
        para <- paste("Log precision for",shrinkpara)
        nch <- nchar(para)
        pxbeta <- lapply(fitall, function(ex2) {
                                mhs <- ex2$internal.marginals.hyperpar
                                if(!is.null(mhs)){
                                nms <- names(mhs) 
                                wh <- which(nms == para)
                                whl <- length(wh)
                                if(whl==0) {
                                nms2 <- sapply(nms,substr,1,nch)
                                wh <- which(nms2==para)
                                }
                                mhswh <- mhs[wh]
                                } else {
                                mhswh <- NULL
                                }
                                return(mhswh)
                                })
    } 
}
rm(fitall,fitall0)
gc()
#whichna <- which(is.na(unlist(lapply(pxbeta,function(pxb) if(length(pxb)==0) NA else pxb[[1]][1]))))

repNA <- function(x) {if(is.na(x)) return(-10^10) else return(x)}  
     
pmt <-proc.time()
if(!only0){
ntag <- length(pxbeta)
funtag <- function(tag){
#if(is.wholenumber(tag/1)) print(tag)
if(is.null(pxbeta[[tag]])) return(NULL) else {
#pbetax <- lapply(as.list(1:20),function(tag){
#tag<-1
postbetanon0all <- list()
postbeta0all <- c()
integralall <- c()
logliksum <- 0
nm <- names(pxbeta[[tag]])
ml0 <- mlik0[tag]
ml1 <- mlik[tag]

ml0 <- repNA(ml0)
ml1 <- repNA(ml1)

maxlik <- max(c(ml0,ml1))

if(!is.null(shrinklc)) whichlc <- which(!(nm %in% shrinkpara)) else whichlc <- NULL
    for(i in 1:length(pxbeta[[tag]])){
    #i<-1        
    pxbetatag <- pxbeta[[tag]][[i]]
    if(!is.null(pxbetatag)) {
        integral <- p0est*exp(ml0-maxlik) + (1-p0est)*exp(ml1-maxlik)
        postbeta0 <- p0est*exp(ml0-maxlik)/integral
        postbetanon0 <- pxbetatag
        postbeta0all <- c(postbeta0all,postbeta0)
        colnames(postbetanon0) <- c("x","y")
    } else {
        postbetanon0 <- NULL
        integral <- 1
    }
    postbetanon0all <- c(postbetanon0all,list(postbetanon0))        
    logliksum<-logliksum+log(integral) + maxlik
    integralall <- c(integralall,integral)  
  
}
names(postbetanon0all) <- nm
postlist<-list(postbetanon0=postbetanon0all,postbeta0=postbeta0all,loglik=logliksum,integralall=integralall)
return(postlist)
}} } else {
ntag <- length(mlik0)
funtag <- function(tag){
#if(is.wholenumber(tag/1)) print(tag)
#pbetax <- lapply(as.list(1:20),function(tag){
#tag<-1
postbeta0all <- c()
integralall <- c()
logliksum <- 0
ml0 <- mlik0[tag]
ml1 <- mlik[tag]

ml0 <- repNA(ml0)
ml1 <- repNA(ml1)

maxlik <- max(c(ml0,ml1))

integral <- p0est*exp(ml0-maxlik) + (1-p0est)*exp(ml1-maxlik)
postbeta0 <- p0est*exp(ml0-maxlik)/integral
postbeta0all <- c(postbeta0all,postbeta0)
logliksum<-logliksum+log(integral) + maxlik
integralall <- c(integralall,integral)  
postlist<-list(postbetanon0=NULL,postbeta0=postbeta0all,loglik=logliksum,integralall=integralall)
return(postlist)
}
}

if(ncpus==1) pbetax <- lapply(as.list(1:ntag),funtag) else {
    mysfExport(forceexport=c("pxbeta"))
    pbetax <- sfLapply(as.list(1:ntag),funtag)
    sfRemoveAll()
    sfStop()
    }
return(pbetax)
}
}
