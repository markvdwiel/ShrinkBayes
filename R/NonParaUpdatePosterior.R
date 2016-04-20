NonParaUpdatePosterior <-
function(fitall, updateoutput, fitall0=NULL, ncpus=2){
#modus: "fixed", "random", "logdisp" 
#shrinkpara: provide when modus is not "logdisp" or "logitp0"
#ndrawtot: total number of draws from the posteriors used to update the prior
#fitall: object fitted using non-shrunken continuous prior
#mufit: mean from fit, usually zero
#priornew: update all posteriors when an estmate of the prior is known. 
#precfitall: scalar or vector. In the latter case different precisions were used for different tags for the initial fitting. 
#gammainit = c(0.001,0.001). Initial values for log-gamma prior; only relevant for modus="random"
#maxsupport: maximum of the support of the prior and posteriors. For numerical stability. -maxsupport is the minimum of the support.
#lcnamederived: name of lincomb components
#fitall <- ShrinkSeq.fitzinb;fitall0<-ShrinkSeq.fitzinb0;fitall0[[1]]<- fitall0[[1]][c(1560,1671)];fitall[[1]]<- fitall[[1]][c(1560,1671)];updateoutput <- SSP0b$prior;ncpus=1;
#fitall <- ShrinkSeq.fitzinb2 ; 

fitall <- fitall[[1]]
if(!is.null(fitall0)){
    fitall0 <- fitall0[[1]]
    zerofit<-TRUE
    } else {
    zerofit <- FALSE
    }
priornew <- updateoutput$priornew
ip <- updateoutput$inputpar
p0  <- updateoutput$p0est
modus       <-ip$modus     
shrinkpara  <-ip$shrinkpara
shrinklc    <-ip$shrinklc  
mufit       <-ip$mufit    
precfit <- ip$precfit
precfitlc <- ip$precfitlc
mufitlc <- ip$mufitlc
gammafit   <-ip$gammafit    
maxsupport  <-ip$maxsupport
pointmass <- ip$pointmass 

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

if(zerofit){
mlik <- mliks(fitall)
mlik0 <- mliks(fitall0)
repNA <- function(x) {if(is.na(x)) return(-10^10) else return(x)}
mlik <- sapply(mlik,repNA)
mlik0 <- sapply(mlik0,repNA)
}

rm(fitall)
if(zerofit) rm(fitall0)
gc()

#whichna <- which(is.na(unlist(lapply(pxbeta,function(pxb) if(length(pxb)==0) NA else pxb[[1]][1]))))

#apply spline to save time
print("Started interpolation")
pmt <- proc.time()
pxbeta <- lapply(1:length(pxbeta), function(i){
        postdisti <- pxbeta[[i]]
     if(is.null(postdisti)) return(NULL) else 
        {
        lapply(postdisti,function(poster){
        if(nrow(poster)<200){
        whpost <- which(abs(poster[,1])<=maxsupport)
        poster <- poster[whpost,]
        prspline1 <- try(myinla.smarginal(poster,len=300))
        if(class(prspline1)=="try-error") return(NULL) else
        {
        prspline2 <- cbind(prspline1$x,prspline1$y)
        return(prspline2)
        }
        } else {
            whpost <- which(abs(poster[,1])<=maxsupport)
            if(length(whpost)<=10) return(NULL) else { #almost no support within required domain
            poster <- poster[whpost,]
            return(poster)
            }
        }
        })
        }})
proc.time()-pmt

if(!zerofit){
pxbeta_eq0 <- lapply(pxbeta, function(postdisti)
    {
        if(!is.null(postdisti)) 
        {
        lapply(postdisti,function(poster)
            {
            if(!is.null(poster))
            {
            pr0 <- myinla.dmarginal(pointmass,poster)
            return(pr0)
            } else {
            return(NULL)
            }})
        } else {NULL}
    })
} 
    
newdens1 <- priornew
nx <- newdens1[,1]
ny <- newdens1[,2]
minx <- max(-maxsupport,nx[1])
maxx <- min(maxsupport,nx[length(nx)])
ndx <- function(ex,nx,ny) {if(ex <= minx) return(0) else 
    {
    if(ex >= maxx) {return(0)} else
    {
    wh <- min(which(nx>=ex))-1
    res <- ny[wh] + ((ny[wh+1]-ny[wh])/(nx[wh+1]-nx[wh]+10^(-10)))*(ex-nx[wh])
    return(res)
    }
    }
}
  
ntag <- length(pxbeta)     
pmt <-proc.time()
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
if(!is.null(shrinklc)) whichlc <- which(!(nm %in% shrinkpara)) else whichlc <- NULL
    for(i in 1:length(pxbeta[[tag]])){
    #i<-1
        if(!is.null(shrinklc) & is.element(i,whichlc)){
        precfiti <- precfitlc 
        mufiti <- mufitlc
        }
        else {
        precfiti <- precfit #divide precision by two for linear combination 
        mufiti <- mufit
        }
        f0init <- dnorm(pointmass,mean=mufiti,sd=1/sqrt(precfiti))
          
    pxbetatag <- pxbeta[[tag]][[i]]
    if(!zerofit) pxbeta_eq0tag <- pxbeta_eq0[[tag]][[i]]
    if(!is.null(pxbetatag)) {
    #whsel0 <- which(abs(pxbetatag[,1])<2)
        support <- pxbetatag[,1]
        
        priornon0 <- sapply(support,ndx,nx=nx,ny=ny)
        if(!zerofit) pxbeta_eq0tagweight <- (pxbeta_eq0tag/f0init)*(p0) else pxbeta_eq0tagweight <- exp(mlik0[tag]-mlik[tag])*p0
        
        #whsel <- which(abs(newdens1[,1])<2)
        supportnd <- newdens1[,1]
        if(modus=="fixed" | modus=="logdisp") finit <- dnorm(supportnd,mean=mufiti,sd=1/sqrt(precfiti)) else
        finit <- dlgamma(supportnd, k=gammafit[1], scale=1/gammafit[2])         
        newdens_sc <- cbind(supportnd,(1-p0)*newdens1[,2]/finit) 
        integral <- pxbeta_eq0tagweight + myinla.expectation_emp(newdens_sc,pxbetatag) 
        if(abs(log(integral)>3)) print(paste("Large integral detected for tag",tag,". May indicate numerical instability."))     
        if(modus=="fixed" | modus=="logdisp") finitsup <- dnorm(support,mean=mufiti,sd=1/sqrt(precfiti)) else
            finitsup <- dlgamma(support, k=gammafit[1], scale=1/gammafit[2])         
        pxbetatagweight <- pxbetatag[,2]*priornon0/finitsup
        #integral2 <- myinla.expectation_emp(cbind(support,priornon0/finitsup),pxbetatag)
        #print(c(integral,integral2))
        postbeta0 <- pxbeta_eq0tagweight/integral
        postbetanon0 <- pxbetatagweight/integral
        postbetanon0 <- cbind(support,postbetanon0)
        postbeta0all <- c(postbeta0all,postbeta0)
        colnames(postbetanon0) <- c("x","y")
    } else {
        postbetanon0 <- NULL
        integral <- 1
    }
    postbetanon0all <- c(postbetanon0all,list(postbetanon0))        
    logliksum<-logliksum+log(integral)
    integralall <- c(integralall,integral)  
  
}
names(postbetanon0all) <- nm
postlist<-list(postbetanon0=postbetanon0all,postbeta0=postbeta0all,loglik=logliksum)
return(postlist)
}
}

if(ncpus==1) pbetax <- lapply(as.list(1:ntag),funtag) else {
    mysfExport(forceexport=c("pxbeta","pxbeta_eq0"),forceexclude=c("fitall","fitall0"))
    pbetax <- sfLapply(as.list(1:ntag),funtag)
    sfRemoveAll()
    sfStop()
    }
logliks <- unlist(lapply(pbetax,function(x) x[[3]]))
lll <- length(logliks)
print(length(logliks[logliks==-Inf]))
logliks <- logliks[logliks!=-Inf]
totalloglik <- mean(logliks)*lll
print(paste("Total log-lik:",totalloglik))
return(pbetax)
} #END FUNCTION
