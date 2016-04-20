NonParaUpdatePrior <-
function(fitall,fitall0=NULL, modus="fixed", shrinkpara=NULL,shrinklc=NULL, lincombs=NULL, includeP0 = TRUE, unimodal=TRUE, 
logconcave=FALSE, symmetric=FALSE, allow2modes=TRUE, maxiter=10, tol=0.005, ndrawtot = 10000, ntotals=c(500,2000,5000,10000), priornew=NULL, 
gaussinit0 = NULL, gammainit0 = NULL, p0init=0.8, p0lower = 0.5, pointmass = 0, maxsupport=6, plotdens=TRUE, ncpus=2){
#
#fitall=fitg;fitall0 <- fitg0;p0init<-0.8;shrinklc = NULL; tol=0.01;maxiter=2;ntotals=c(2000);modus="fixed"; shrinkpara="group"; ncpus=6;
#shrinklc=FALSE;ndrawtot = 10000;priornew=NULL;p0init=0.7; p0lower = 0.5; pointmass = 0;spikeslab=FALSE;
#gaussinit0=NULL; gammainit0 = NULL; maxsupport=5; plotdens=FALSE;unimodal=TRUE;includeP0 = TRUE;logconcave=TRUE;symmetric=FALSE; allow2modes=TRUE;
#maxsupport <- 5
if(is.null(shrinkpara) & is.null(shrinklc)) {
print("PLEASE SPECIFY EITHER OF THE ARGUMENTS shrinkpara OR shrinklc")
return(NULL)
}

if(!includeP0) p0init <- 0
paramprior <- fitall[[2]]
fitall <- fitall[[1]]
mufit<-paramprior$mufixed
precfit<-paramprior$precfixed

maxsupport <- max(maxsupport,6/precfit)

#muaddinit<-paramprior$muaddfixed
#precaddfit<-paramprior$precaddfixed
if(!is.null(fitall0)){
    fitall0 <- fitall0[[1]]
    zerofit<-TRUE
    } else {
    zerofit <- FALSE
    }

if(!is.null(shrinklc)) {
    if(is.null(lincombs)) linc <- fitall[[1]]$.args$lincomb else linc <- lincombs
    if(is.null(linc)) print("Coefficients of linear combination(s) are unknown. Please provide the linear combinations.")
    whlc <- which(names(linc) %in% shrinklc) 
    lccoef <- unlist(linc[whlc[1]]) #assume same set of weights for all lc!
    precfitlc <- 1/sum(sapply(1:length(lccoef), function(i){
    coefi <- lccoef[i]
    return(coefi^2/precfit)
    })) 
    mufitlc <- sum(sapply(1:length(lccoef), function(i){
    coefi <- lccoef[i]
    return(coefi*mufit)
    })) 
} else {
mufitlc <- NULL
precfitlc <- NULL
lccoef <- NULL
}

gammainit <- c(paramprior$shaperand,paramprior$raterand)

inputpar = list(modus=modus,shrinkpara=shrinkpara, shrinklc=shrinklc, unimodal=unimodal, 
logconcave=logconcave, symmetric=symmetric,maxiter=maxiter, tol=tol, ndrawtot = ndrawtot, ntotals=ntotals,priorold=priornew,
mufit=mufit,mufitlc=mufitlc, precfit=precfit, precfitlc=precfitlc,gammafit = gammainit,lccoef=lccoef,maxsupport = maxsupport,
gaussinit0 = gaussinit0, gammainit0 = gammainit0, p0init=0.8, pointmass = 0)

if(modus=="fixed") {if(!is.null(shrinkpara)){if(is.factor(try(get(shrinkpara),silent=TRUE))) shrinkpara <- fact2vec(shrinkpara)}}

if(ncpus > 1) {
    sfInit(parallel=TRUE,cpus=ncpus) 
    sfLibrary(VGAM)
    sfLibrary(ShrinkBayes)
    }

if(modus=="fixed" | modus =="logdisp"){
    if(!is.null(gaussinit0) & is.null(priornew)) {
        muinit0 <- gaussinit0[1];precinit0 <- gaussinit0[2]
        sup <- muinit0 + seq(-maxsupport,maxsupport,by=0.05)
        priornew <- cbind(sup,dnorm(sup,mean=muinit0,sd=sqrt(1/(precinit0))))
        }
        
    if(is.null(gaussinit0) & is.null(priornew)) {
        sup <- mufit + seq(-maxsupport,maxsupport,by=0.05)
        priornew <- cbind(sup,dnorm(sup,mean=mufit,sd=sqrt(1/(precfit))))
        }
}

if(modus=="random"){
    if(!is.null(gammainit0) & is.null(priornew)) {
        sup <- seq(-maxsupport,maxsupport,by=0.05)
        priornew <- cbind(sup,dlgamma(sup,k=gammainit0[1], scale=1/gammainit0[2]))
        }
        
    if(is.null(gammainit0) & is.null(priornew)) {
        sup <- seq(-maxsupport,maxsupport,by=0.05)
        priornew <- cbind(sup,dlgamma(sup,k=gammainit[1], scale=1/gammainit[2]))
        }
}


p0 <- p0init    
lntotal <- length(ntotals)
densall <- list()
ksmaxall <- c()
totalloglikall <- c()
relchangeall <- c()
sdall <-c()
p0all <- c()
nall <- c()


d <- 1
contn <- TRUE
while(d<=lntotal & contn){
#d<-1
ntotal <- ntotals[d]
print(paste("Using (at most)",ntotal,"posteriors")) 



wh <- which.min(unlist(lapply(fitall,is.null)))
callmode <- fitall[[wh]]$call

if(is.null(callmode)) { #input is NOT AN inla output object
pxbeta <- fitall  
} else {
    if(modus=="fixed") { #input is inla output object
pxbeta <- lapply(fitall, function(ex2) {
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

repNA <- function(x) {if(is.na(x)) return(-10^10) else return(x)}
mlikall <- mliks(fitall)
mlikall <- sapply(mlikall,repNA)
if(zerofit){
mlik0all <- mliks(fitall0)
mlik0all <- sapply(mlik0all,repNA)
}




whichna <- which(is.na(unlist(lapply(pxbeta,function(pxb) if(length(pxb)==0) NA else pxb[[1]][1]))))

if(length(whichna)>0){
    pxbeta <- pxbeta[-whichna]
    mlikall <- mlikall[-whichna]
    if(zerofit){
        mlik0all <- mlik0all[-whichna]
    }
}



if(!is.null(priornew)) {
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
    }

ntag0<- length(pxbeta) 
ntag <- ntag0
ncomp <- length(pxbeta[[1]])
ntot <- ntag*ncomp

arr <- matrix(0, ncol=ncomp,nrow=ntag0)
if(ntotal < (ntag0*ncomp)) {
ntot <- ntotal
set.seed(74574)
wh <- sample(1:(ncomp*ntag0),ntotal)
} else {
    wh <- 1:(ntag0*ncomp)
    contn <- FALSE  #end iteration when ntotal exceeds the total number of available posteriors
}

which1 <- sapply(wh, function(whi) {
#whi <- 1904
arra <- arrayInd(whi,c(ntag0,ncomp))
pxbetaij <- pxbeta[[arra[1]]][[arra[2]]]
if(!is.null(pxbetaij)){
    nr <- nrow(pxbetaij)
    lpost <- length(which(abs(pxbetaij[,1])<=maxsupport))
    if (lpost/nr <= 0.1) return(0) else return(1)
} else return(0)
})
arr[wh] <- which1
whzero <- which(apply(arr,1,sum)==0)
if(length(whzero)>0) {
    pxbeta <- pxbeta[-whzero]
    if(!is.null(fitall0)){
        mlik <- mlikall[-whzero]
        mlik0 <- mlik0all[-whzero]
    }
    arr <- arr[-whzero,,drop=FALSE]
    } else {
    if(!is.null(fitall0)){
        mlik <- mlikall
        mlik0 <- mlik0all
        }
    }
ntag <- nrow(arr)
for(i in 1:ntag){
pxbeta[[i]] <- pxbeta[[i]][which(arr[i,]==1)]
}



#apply spline to save time
pmt <- proc.time()
print("Started interpolation")
pxbeta <- lapply(pxbeta, function(postdisti){lapply(postdisti,function(poster){
    if(nrow(poster)<200){
    whpost <- which(abs(poster[,1])<=maxsupport)
    poster <- poster[whpost,]
    prspline1 <- try(myinla.smarginal(poster,len=300))
    if(class(prspline1)=="try-error") return(poster) else {
    prspline2 <- cbind(prspline1$x,prspline1$y)
    return(prspline2)
    }
    } else { #note rescaling not necessary because all distributional properties are computed using spline and internal standardization
        whpost <- which(abs(poster[,1])<=maxsupport)
        poster <- poster[whpost,]
        return(poster)   
        }
    })})
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
    })} else {
pxbeta_eq0 <-NULL
}
        
if(is.null(p0init)) {
p0init<-0.8
print(paste("No initial p0 determined. Set to",p0init))
}

    
iter<-1
ksmax <- 1
ndraw <- ceiling(ndrawtot/ntot * (1/max(0.01,(1-p0))))
totalloglikprev <- 0  

integralall <- c()
likinc <- TRUE
while(iter <= maxiter & ksmax >= tol & likinc){
pmt <-proc.time()
    print(paste("iteration:",iter))
    funtag <- function(tag){
    #if(is.wholenumber(tag/50)) print(tag)
    #tag <- 1
    #print(tag)
    postbetanon0all <- list()
    postbeta0all <- c()
    logliksum <- 0
    if(!is.null(shrinklc)) whichlc <- which(!(names(pxbeta[[tag]]) %in% shrinkpara)) else whichlc <- NULL
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
            
        support <- pxbetatag[,1]
        if (iter==1 & is.null(priornew)){
            pxbetatagweight <- pxbetatag[,2]
            integral <-1-p0
            if(!zerofit) pxbeta_eq0tagweight <- (pxbeta_eq0tag/f0init)*(p0) else pxbeta_eq0tagweight <- exp(mlik0[tag]-mlik[tag])*p0
            integral <- integral + pxbeta_eq0tagweight
            } else {
            priornon0 <- sapply(support,ndx,nx=nx,ny=ny)
            if(!zerofit) pxbeta_eq0tagweight <- (pxbeta_eq0tag/f0init)*(p0) else pxbeta_eq0tagweight <- exp(mlik0[tag]-mlik[tag])*p0
            supportnd <- newdens1[,1]
            if(modus=="fixed" | modus=="logdisp") finit <- dnorm(supportnd,mean=mufiti,sd=1/sqrt(precfiti)) else
            finit <- dlgamma(supportnd, k=gammainit[1], scale=1/gammainit[2])         
            newdens_sc <- cbind(supportnd,(1-p0)*newdens1[,2]/finit) 
            
            integral <- pxbeta_eq0tagweight + myinla.expectation_emp(newdens_sc,pxbetatag)
          }
        if(modus=="fixed" | modus=="logdisp") finitsup <- dnorm(support,mean=mufiti,sd=1/sqrt(precfiti)) else
            finitsup <- dlgamma(support, k=gammainit[1], scale=1/gammainit[2])         
        if(iter==1 & is.null(priornew)) pxbetatagweight <- pxbetatag[,2] else pxbetatagweight <- pxbetatag[,2]*priornon0/finitsup
        postbeta0 <- pxbeta_eq0tagweight/integral
        postbetanon0 <- pxbetatagweight/integral
        postbetanon0 <- cbind(support,postbetanon0)
        colnames(postbetanon0) <- c("x","y")
        postbeta0all <- c(postbeta0all,postbeta0)
        postbetanon0all <- c(postbetanon0all,list(postbetanon0))        
        logliksum<-logliksum+log(integral)
        integralall <- c(integralall,integral)     
    }   
      
    samples<-lapply(1:length(postbetanon0all),function(i){
        post1 <- postbetanon0all[[i]]
        samps <- samplepost(post1,ndraw)
        samps <- samps[abs(samps)<=maxsupport]
        return(samps)
        })
        
    finsamps <- sapply(1:length(postbeta0all),function(j) {
        post0 <- postbeta0all[j]
        samx <- samples[[j]]
        if(length(samx)>0){
        onezerosamx <- sapply(samx,function(x){
        if(post0 > 0 & post0<1){ 
            onezero <- sample(c(x,0),size=1,prob=c(1-post0,post0))
        } else { 
            if(post0==0){
                onezero <- x
            } else {
                onezero <- 0
            }
        }
        return(onezero)
        })
        return(onezerosamx)
        } else return(list(c()))
        })
        
    postlist=list(samples=finsamps,postbeta0=postbeta0all,loglik=logliksum,integralall=integralall)
    return(postlist)
    }
    if(ncpus==1) pbetax <- lapply(as.list(1:ntag),funtag) else {
    if(iter==1) {
        print("Exporting to slaves") 
        mysfExport(forceexport=c("pxbeta","pxbeta_eq0"),forceexclude=c("fitall","fitall0"))
        print("Finished exporting")
        } else {
        sfExport("ndx","newdens1","p0","priornew","iter")
        }
    pbetax <- sfLapply(as.list(1:ntag),funtag)
    }
    
    integ <- unlist(lapply(pbetax,function(x) x[[4]]))
    logliks <- unlist(lapply(pbetax,function(x) x[[3]]))
    lll <- length(logliks)
    print(length(logliks[logliks==-Inf]))
    logliks <- logliks[logliks!=-Inf]
    totalloglik <- mean(logliks)*lll
    relchange <- abs((totalloglik - totalloglikprev)/totalloglik)
    relchangenoabs <- (totalloglik - totalloglikprev)/totalloglik
    totalloglikprev <- totalloglik
    
    meanp0<-mean(unlist(lapply(pbetax,function(x) x[[2]]))) 
    postfixed <- unlist(lapply(pbetax,function(x) x[[1]])) 
        
    
  
    print(paste("Prop. NA:",length(postfixed[is.na(postfixed)])/length(postfixed)))
    postfixed <- postfixed[!is.na(postfixed)]
    postfixed <- postfixed[postfixed!=0]
    if(symmetric) postfixed <- c(postfixed,2*median(postfixed)-postfixed)
    
    if(allow2modes){
    r <- diff(range(postfixed))
        mygrid <- seq(min(postfixed) - 0.1 * r, max(postfixed) + 0.1 * r, length = 500)

    for(j in 1:2){
        if(j==1) {
        postfixedone <- postfixed[postfixed<0]  
        } else {
        postfixedone <- postfixed[postfixed>0]
        }
        pone <- length(postfixedone)/length(postfixed)
        newdens <- density(postfixedone)
        print("Fitting unimodal density")
        if(iter==1) library(Iso)
        newdens_iso <- ufit(y=newdens$y,x=newdens$x)
        newdens$x <- newdens_iso$x
        newdens$y <- newdens_iso$y
        print("Fitting done")
        if(j==1) {newy <- pone*inla.dmarginal(mygrid,cbind(newdens$x,newdens$y))}
        if(j==2) {newy <- pone*inla.dmarginal(mygrid,cbind(newdens$x,newdens$y)) + newy}
        }
        newdens <- list(x=mygrid,y=newy)
    } else {
    newdens <- density(postfixed)
    #newdensiso <- newdens
    if(unimodal & !logconcave){
    print("Fitting unimodal density")
    if(iter==1) library(Iso)
    newdens_iso <- ufit(y=newdens$y,x=newdens$x)
    newdens$x <- newdens_iso$x
    newdens$y <- newdens_iso$y
    print("Fitting done")
    }
    
    if(logconcave){
    if(iter==1) library(logcondens)
    print("Fitting log-concave density")
    res <- logConDens(postfixed, smoothed = TRUE, print = FALSE)
    newdens$x <- res$xs
    newdens$y <- res$f.smoothed
    print("Fitting done")
    }  
    }

    
    
    if(plotdens) {
    plot(newdens)
    }
     
    nx <- newdens$x
    ny <- newdens$y
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
    
    sdest <- sd(postfixed)        
    newdens1 <- cbind(nx,ny)
    if(includeP0) p0 <- max(p0lower,meanp0) else p0 <- 0
    if(iter==1) postfixedprev <- rnorm(ndrawtot,mufit,1/sqrt(precfit)) 
     
    ksmax <- ks.test(postfixedprev,postfixed, alternative = "two.sided")$statistic
    postfixedprev <- postfixed 
    densall <- c(densall,list(newdens1))
    ksmaxall <- c(ksmaxall,ksmax)
    sdall <- c(sdall,sdest)
    p0all <- c(p0all,p0)
    nall <- c(nall,ntotal)
    totalloglikall <- c(totalloglikall,totalloglik)
    if((iter > 2) & relchangenoabs<0) likinc <- FALSE
    print(likinc)
    iter <- iter+1
    relchangeall <- c(relchangeall,relchange) 
    conv <- c(ksmax,totalloglik,relchange,sdest)
    names(conv) <- c("ksmax","totalloglik","relchange","sd")
    print(conv)
    print(paste("Current p0 estimate:",p0))
    proc.time()-pmt
} #END INNER WHILE LOOP

if(ncpus>1) {
    sfRemoveAll()
    }
whncur <- which(nall==ntotal)
whmax <- which.max(totalloglikall[whncur])
if(whmax>1) {
    priornew <- densall[whncur][[whmax-1]] #previous parameter configuration corresponds to maximum marg lik
    p0 <- p0all[whncur][whmax-1]
    }
d <- d+1
} #END OUTER WHILE LOOP

if(!is.null(gaussinit0) & is.null(priornew)) sd1 <- 1/sqrt(gaussinit0[2]) else sd1 <- NA
densall <- c(list("Initial prior"),densall);ksmaxall <- c(NA,ksmaxall);sdall <- c(sd1,sdall); p0all <- c(p0init,p0all);
totalloglikall<-c(totalloglikall,NA);nallforlik<-c(nall,NA)

nmax <- max(nallforlik,na.rm=T)
whsel <- which(nallforlik==nmax)
#print(whsel)
priors <- densall[whsel]
tll <- totalloglikall[whsel]
wh <- which.max(tll)
#print(wh)
prior <- priors[[wh]] #FINAL PRIOR IS THE ONE WITH THE LARGEST MARGINAL LIKELIHOOD
p0est <- p0all[whsel][wh]
toret <- list(priornew=prior, p0est = p0est, inputpar = inputpar, densall=densall,p0all=p0all,ksmaxall=ksmaxall,sdall = sdall, totalloglikall=totalloglikall,nallforlik=nallforlik,notfit=whichna)
return(toret)
} #END FUNCTION
