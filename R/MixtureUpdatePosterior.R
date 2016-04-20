MixtureUpdatePosterior <-
function(fitall, updateoutput, fitall0=NULL, ncpus=1){
#fitall <- fitg;updateoutput<-prior;fitall0=NULL;ncpus=8
#updateoutput<- mixtprior
fitall <- fitall[[1]]
if(!is.null(fitall0)){
    fitall0 <- fitall0[[1]]
    zerofit<-TRUE
    } else {
    zerofit <- FALSE
    }

best <- updateoutput$best
nam <- names(best)
p0 <- best[which(nam=="p0")]
stdev <- best[which(nam=="stdev")]

ip <- updateoutput$inputpar
modus <-ip$modus  
if(modus=="mixt") {
    pminus <- best[which(nam=="pmin")] 
    mu <- best[which(nam=="mu")]
    }  
shrinkpara  <-ip$shrinkpara
shrinklc    <-ip$shrinklc  
mufit      <-ip$mufit 
precfit  <-ip$precfit
precfitlc  <-ip$precfitlc   
mufitlc <- ip$mufitlc 
maxsupport  <-ip$maxsupport
pointmass <- ip$pointmass


if(!is.null(shrinkpara)){if(is.factor(try(get(shrinkpara),silent=TRUE))) shrinkpara <- fact2vec(shrinkpara)}




callmode <- fitall[[1]]$call

if(!is.null(callmode)) { #input is inla output object
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
} else { #input is list of posteriors
pxbeta <- fitall
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
#
#
#if(length(whichna)>0){
#    #postdist<-postdist[-whichna]
#    pxbeta <- pxbeta[-whichna]
#    mlik <- mlik[-whichna]
#    mlik0 <- mlik0[-whichna]
#}

if(ncpus > 1) {
    sfInit(parallel=TRUE,cpus=ncpus) 
    sfLibrary(ShrinkBayes)
    }



#apply spline to save time
pmt <- proc.time()
print("Started interpolation")
pxbeta <- lapply(1:length(pxbeta), function(i){
        postdisti <- pxbeta[[i]]
     if(is.null(postdisti)) return(NULL) else 
        {
        lapply(postdisti,function(poster){
        if(nrow(poster)<200){
        whpost <- which(abs(poster[,1])<=maxsupport)
        poster <- poster[whpost,]
        prspline1 <- try(myinla.smarginal(poster,len=300),silent=T)
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
            pr0 <- try(myinla.dmarginal(pointmass,poster),silent=T)
            if(class(pr0)=="try-error") return(NULL) else return(pr0)
            } else {
            return(NULL)
            }})
        } else {NULL}
    })
}


loglikall <- c()
ntagcur <- length(pxbeta)
pmt <- proc.time()
    tagfun <- function(tag){
    if(is.null(pxbeta[[tag]])) return(NULL) else {
     if(is.wholenumber(tag/500)) print(paste("Tag:",tag))
        #tag<-333
        #print(tag)
        postbetanon0all <- list()
        postbeta0all <- c()  
        logliksum<-0
        nm <- names(pxbeta[[tag]])
        if(!is.null(shrinklc)) whichlc <- which(!(nm %in% shrinkpara)) else whichlc <- NULL
        for(i in 1:length(pxbeta[[tag]])){
        #tag<-1;i<-1
            if(!is.null(shrinklc) & is.element(i,whichlc)){
            precfiti <- precfitlc 
            mufiti <- mufitlc 
            } else {
            precfiti <- precfit #divide precision by two for linear combination 
            mufiti <- mufit
            }
            f0init <- dnorm(0,mean=mufiti,sd=1/sqrt(precfiti))
            pxbetatag <- pxbeta[[tag]][[i]]
            if(!zerofit) {if(is.null(pxbeta_eq0[[tag]][[i]])) pxbeta_eq0tag <- 10^6 else pxbeta_eq0tag <- pxbeta_eq0[[tag]][[i]]}
            fxinit <- function(x) dnorm(x,mean=mufiti,sd=1/sqrt(precfiti)) 
            support <- pxbetatag[,1]
            finit <-dnorm(support,mean=mufiti,sd=1/sqrt(precfiti))
            if(modus=="gauss") {
                if(!zerofit) integral0 <- (pxbeta_eq0tag/f0init)*(p0) else integral0 <- exp(mlik0[tag]-mlik[tag])*p0
      
                integralnon0 <- myinla.expectation(function(x) (1-p0)*dnorm(x,mean=0,sd=stdev)/fxinit(x),marginal=pxbetatag)
                integral <- integral0 + integralnon0
                pxbetatagweight <- (pxbetatag[,2]/finit)* (1-p0)*dnorm(support,mean=0,sd=stdev)
            }
            if(modus=="laplace") {
                sc <- stdev/sqrt(2)
                if(!zerofit) integral0 <- (pxbeta_eq0tag/f0init)*(p0) else integral0 <- exp(mlik0[tag]-mlik[tag])*p0
                integralnon0 <- myinla.expectation(function(x) (1-p0)*dlaplace(x,location=0,scale=sc)/fxinit(x),marginal=pxbetatag)
                integral <- integral0 + integralnon0 
                pxbetatagweight <- (pxbetatag[,2]/finit)* (1-p0)*dlaplace(support,location=0,scale=sc)
            }
            #if(modus=="gammamixt") {
#                ra <- mu/stdev^2
#                sh <- mu^2/stdev^2
#                if(!zerofit) integral0 <- (pxbeta_eq0tag/f0init)*(p0) else integral0 <- exp(mlik0[tag]-mlik[tag])*p0
#                integralnon0 <- myinla.expectation(function(x) (1-p0)*((1-pwide)*(1/2*dgamma(-x,shape=sh,rate=ra) + 1/2*dgamma(x,shape=sh,rate=ra)) + 
#                pwide*dnorm(x,mean=0,sd=sdwide))/fxinit(x),marginal=pxbetatag) 
#                integral <- integral0 + integralnon0 
#                pxbetatagweight <- (pxbetatag[,2]/finit)* (1-p0)*((1-pwide)*(1/2*dgamma(-support,shape=sh,rate=ra) + 1/2*dgamma(support,shape=sh,rate=ra)) + 
#                pwide*dnorm(support,mean=0,sd=sdwide))
#            }
            if(modus=="mixt") {
                if(!zerofit) integral0 <- (pxbeta_eq0tag/f0init)*(p0) else integral0 <- exp(mlik0[tag]-mlik[tag])*p0
                integralnon0 <- myinla.expectation(function(x) (1-p0)*(pminus*dnorm(x,mean=-mu,sd=stdev) + (1-pminus)*dnorm(x,mean=mu,sd=stdev))/fxinit(x),marginal=pxbetatag)  
                integral <- integral0 + integralnon0 
                pxbetatagweight <- (pxbetatag[,2]/finit)* (1-p0)*(pminus*dnorm(support,mean=-mu,sd=stdev) + (1-pminus)*dnorm(support,mean=mu,sd=stdev))
            }
            postbetanon0 <- pxbetatagweight/integral
            if(p0!=0) postbeta0 <- integral0/integral else postbeta0 <- 0
            postbetanon0 <- cbind(support,postbetanon0)
            colnames(postbetanon0) <- c("x","y")
            postbetanon0all <- c(postbetanon0all,list(postbetanon0))       
            postbeta0all <- c(postbeta0all,postbeta0)  
            logliksum<-logliksum+log(integral) 
        }
        names(postbetanon0all) <- names(postbeta0all) <- nm
        postlist=list(postbetanon0=postbetanon0all,postbeta0=postbeta0all,loglik=logliksum)
        return(postlist)   
    }
    }
    
tagfuntry <- function(ind) {
    tr <- try(tagfun(ind),silent=T)
    if(class(tr)=="try-error") return(NULL) else return(tr)
    }
if(ncpus==1) reslogliks <- lapply(1:ntagcur,tagfuntry) else {
mysfExport(forceexport=c("pxbeta","pxbeta_eq0"))
reslogliks <- sfLapply(1:ntagcur,tagfuntry)
}
return(reslogliks)
} #END FUNCTION
