MixtureUpdatePrior <-
function(fitall,fitall0=NULL, shrinkpara=NULL, modus="mixt", shrinklc=NULL,  lincombs=NULL,ntotal = 10000, maxsupport=6,pointmass=0,
pminvec = c(0,0.25,0.5,0.75,1), p0vec = c(0.5,0.7,0.9,0.95,0.99,0.999,1), 
meanvec = c(0.1, 0.3, 0.5, 0.75,1.5),sdvec=c(0.2,0.5,0.8,1.5,3),meansdauto=TRUE,
ncpus=2,refinegrid=TRUE, symmetric=FALSE){

# # 
 fitall <- fitg; fitall0=fitg0;shrinkpara="groupsIG"; modus="gauss"; ntotal = 10000;
 maxsupport=6;pminvec=c(0.1,0.3,0.5);pointmass=0;lincombs <- NULL;shrinklc=NULL;
 p0vec = c(0.5,0.8,0.9,0.99,1); meanvec = c(0.1, 0.25, 0.4, 0.5, 0.75);sdvec=c(0.2,0.75,1.5);ncpus=6;refinegrid=TRUE;
 meansdauto=TRUE; symmetric=FALSE
if(is.null(shrinkpara) & is.null(shrinklc)) {
print("PLEASE SPECIFY EITHER OF THE ARGUMENTS shrinkpara OR shrinklc")
return(NULL)
}
paramprior <- fitall[[2]]
fitall <- fitall[[1]]
muinit=paramprior$mufixed
precfit=paramprior$precfixed

if(symmetric) pminvec <- c(0.5)

maxsupport <- max(maxsupport,6/precfit)

if(!is.null(fitall0)){
    fitall0 <- fitall0[[1]]
    zerofit<-TRUE
    } else {
    zerofit <- FALSE
    }

#muaddinit<-paramprior$muaddfixed
#precaddinit<-paramprior$precaddfixed
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
    return(coefi*muinit)
    })) 
} else {
mufitlc <- NULL
precfitlc <- NULL
lccoef <- NULL
}

inputpar <- list(modus, shrinkpara, shrinklc,lincombs, ntotal, maxsupport,pointmass,pminvec,
p0vec, meanvec,sdvec,muinit,mufitlc,precfit,precfitlc)
names(inputpar) <- c("modus","shrinkpara","shrinklc","lincombs","ntotal","maxsupport","pointmass","pminvec","p0vec","meanvec",
"sdvec","mufit","mufitlc","precfit","precfitlc")

if(!is.null(shrinkpara)){if(is.factor(try(get(shrinkpara),silent=TRUE))) shrinkpara <- fact2vec(shrinkpara)}


callmode <- fitall[[1]]$call

if(!is.null(callmode)) { #input is inla output object
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
} else { #input is list of posteriors
pxbeta <- fitall
}

repNA <- function(x) {if(is.na(x)) return(-10^10) else return(x)}
mlikall <- mliks(fitall)
mlikall <- sapply(mlikall,repNA)
if(zerofit){
mlik0all <- mliks(fitall0)
mlik0all <- sapply(mlik0all,repNA)
#if(is.null(p0init)) p0init <- mean(exp(mlik0all-mlikall)/(1+exp(mlik0all-mlikall)))
#print(paste("Initial p0",p0init))
}

whichna <- which(is.na(unlist(lapply(pxbeta,function(pxb) if(length(pxb)==0) NA else pxb[[1]][1]))))


if(length(whichna)>0){
    #postdist<-postdist[-whichna]
    pxbeta <- pxbeta[-whichna]
    mlikall <- mlikall[-whichna]
    if(zerofit){
        mlik0all <- mlik0all[-whichna]
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
nr <- nrow(pxbetaij)
lpost <- length(which(abs(pxbetaij[,1])<=maxsupport))
if (lpost/nr <= 0.1) return(0) else return(1)
})
arr[wh] <- which1
whzero <- which(apply(arr,1,sum)==0)
if(length(whzero)>0) {
    pxbeta <- pxbeta[-whzero]
    arr <- arr[-whzero,,drop=FALSE]
    mlikall <- mlikall[-whzero]
    if(zerofit) mlik0all <- mlik0all[-whzero]
    }
ntag <- nrow(arr)
for(i in 1:ntag){
pxbeta[[i]] <- pxbeta[[i]][which(arr[i,]==1)]
}

if(length(whzero)>0) tagsinv <- (1:ntag0)[-whzero] else tagsinv <- (1:ntag0)

#apply spline to save time
pmt <- proc.time()
print("Started interpolation")
pxbeta <- lapply(pxbeta, function(postdisti){lapply(postdisti,function(poster){
    if(nrow(poster)<200){
    whpost <- which(abs(poster[,1])<=maxsupport)
    poster <- poster[whpost,]
    prspline1 <- myinla.smarginal(poster,len=300)
    prspline2 <- cbind(prspline1$x,prspline1$y)
    return(prspline2)
    } else {
        whpost <- which(abs(poster[,1])<=maxsupport)
        poster <- poster[whpost,]
        return(poster)
    }})})
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

if(!zerofit) {
if(is.null(shrinklc)) f0init <- dnorm(0,mean=muinit,sd=1/sqrt(precfit)) else f0init <- dnorm(0,mean=mufitlc,sd=1/sqrt(precfitlc))
mlik0all <- log(unlist(pxbeta_eq0)/f0init) + mlikall
}

#for numerical stability; rationale is that the null-model can not have a log-marg lik > log-mark lik (full model) + 5; for numerical stability
difffun <- function(diff) return(min(5,max(diff,-500)))
mlik0minmlik1 <- sapply(mlik0all-mlikall,difffun)

if(meansdauto){
meanest <- unlist(lapply(pxbeta, function(postdisti)
    {
        if(!is.null(postdisti)) 
        {
        lapply(postdisti,function(poster)
            {
            if(!is.null(poster))
            {
            #pr0 <- INLA:::inla.expectation(function(x) x,poster)
              pr0 <- myinla.expectation(function(x) x,poster)
            return(pr0)
            } else {
            return(NULL)
            }})
        } else {NULL}
    }))

    
p0init <- min(0.98,mean(exp(mlik0minmlik1)/(1+exp(mlik0minmlik1))))
absmean <- abs(meanest)
whsm <- which(absmean >= quantile(absmean,p0init))
meansel <- meanest[whsm]
if(modus!="mixt") {
    sdest <- mad(meansel)
    print(paste("sdest initial:",sdest))
    sdvec <- c(0.5*sdest,0.75*sdest,sdest,1.5*sdest,2*sdest,3*sdest)
    }
if(modus=="mixt") {
    meanselpos <- meansel[meansel>0]
    meanselneg <- meansel[meansel<0]
    meanest <- median(c(meanselpos,-meanselneg))
    print(paste("mean est initial:",meanest))
    meanvec <- c(0.5*meanest,0.75*meanest,meanest,1.5*meanest,2*meanest)
    sdest <-  mad(c(meanselpos[meanselpos>meanest]-meanest,-meanselpos[meanselpos>meanest]+meanest,
    meanselneg[meanselneg < -meanest]+meanest,-(meanselneg[meanselneg < -meanest]+meanest)))
    print(paste("sd est initial:",sdest))
    sdvec <- c(0.5*sdest,sdest,1.5*sdest,3*sdest)
}
}


if(ncpus > 1) {
    sfInit(parallel=TRUE,cpus=ncpus) 
    sfLibrary(ShrinkBayes)
    }

if(refinegrid) maxrefine <- 2 else maxrefine <- 1
for(refine in 1:maxrefine){

if(modus=="mixt") {
gridd <- c(); for(i0 in 1:length(pminvec)){for (i in 1:length(p0vec)) {
for (j in 1:length(meanvec)) {for (k in 1:length(sdvec)) gridd <- rbind(gridd,c(pminvec[i0],p0vec[i],meanvec[j],sdvec[k]))}}}
colnames(gridd) <- c("pmin","p0","mu","stdev")
} else {
gridd <- cbind(as.vector(apply(matrix(p0vec,ncol=1),1,function(x) rep(x,length(sdvec)))),rep(sdvec,length(p0vec)))
colnames(gridd) <- c("p0","stdev")
}

   

loglikall <- c()
ngrid <- nrow(gridd)
print(paste("nrow grid",ngrid))

print(paste("ntag used",ntag))

ntagmin <- 500;ngridmin <- 25;
ntagreeks <- c()
ngridreeks <- c()
if(ntag<=ntagmin | ngrid<=25){ntagreeks <- c(ntag);ngridreeks <- c(ngrid)} else {
    ntagnew <- 500
    ngridnew <- ngrid
    while(ntagnew<ntag & ngridnew>25){
        ntagreeks <- c(ntagreeks,ntagnew)
        ngridreeks <- c(ngridreeks,max(ngridnew,25))
        ntagnew <- ntagnew*2
        ngridnew <- ceiling(ngridnew/2)
        }
    if(ntag > 1.25*ntagreeks[length(ntagreeks)]) {
        ntagreeks <- c(ntagreeks,ntag)
        ngridreeks <- c(ngridreeks,max(ngridnew,25))
        }
}

for(ng in 1:length(ntagreeks)){
#ng <-1
ntagcur <- ntagreeks[ng]
ngrid<-ngridreeks[ng]
print(ntagcur)
print(ngrid)
pmt <- proc.time()
    tagfun <- function(tag){
     #if(is.wholenumber(tag/50)) print(paste("Tag:",tag))
     #print(paste("Tag:",tag))
        #tag<-9999
        postbetanon0all <- list()
        postbeta0all <- c()  
        logliksum<-0
         
        if(!is.null(shrinklc)) whichlc <- which(!(names(pxbeta[[tag]]) %in% shrinkpara)) else whichlc <- NULL
   
        for(i in 1:length(pxbeta[[tag]])){
        #tag<-3;i<-1
        if(!is.null(shrinklc) & is.element(i,whichlc)){
        precfiti <- precfitlc 
        mufiti <- mufitlc
        } else {
        precfiti <- precfit 
        mufiti <- muinit
        }
            f0init <- dnorm(0,mean=mufiti,sd=1/sqrt(precfiti))
            pxbetatag <- pxbeta[[tag]][[i]]
            if(!zerofit) pxbeta_eq0tag <- pxbeta_eq0[[tag]][[i]]
            fxinit <- function(x) dnorm(x,mean=mufiti,sd=1/sqrt(precfiti)) 
            if(modus=="gauss") {
                fungrid <- function(grrow){
                p0 <- gridd[grrow,1]
                sdzero2<- gridd[grrow,2] 
                #integral0 <- (pxbeta_eq0tag/f0init)*(p0)
                if(!zerofit) integral0 <- (pxbeta_eq0tag/f0init)*(p0) else integral0 <- exp(mlik0minmlik1[tag])*p0
                integralnon0 <- myinla.expectation(function(x) (1-p0)*dnorm(x,mean=0,sd=sdzero2)/fxinit(x),marginal=pxbetatag)
                integral <- integral0 + integralnon0 
                if(log(integral)>3) print(c(grrow,tag,i,log(integral)))
                return(log(integral))
                }
                resall <- sapply(1:ngrid,fungrid) 
            }
            if(modus=="laplace") {
                resall <- sapply(1:ngrid,function(grrow){
                p0 <- gridd[grrow,1]
                sc <- gridd[grrow,2]/sqrt(2)
                if(!zerofit) integral0 <- (pxbeta_eq0tag/f0init)*(p0) else integral0 <- exp(mlik0minmlik1[tag])*p0
                integralnon0 <- myinla.expectation(function(x) (1-p0)*dlaplace(x,location=0,scale=sc)/fxinit(x),marginal=pxbetatag)
                integral <- integral0 + integralnon0 
                return(log(integral))
                })
            }
            #if(modus=="gamma") {
#                resall <- sapply(1:ngrid,function(grrow){ #grrow <- 1
#                p0 <- gridd[grrow,1]
#                mu <- gridd[grrow,2]
#                va <- (gridd[grrow,3])^2
#                ra <- mu/va
#                sh <- mu^2/va
#                if(!zerofit) integral0 <- (pxbeta_eq0tag/f0init)*(p0) else integral0 <- exp(mlik0all[tag]-mlikall[tag])*p0
#                integralnon0 <- myinla.expectation(function(x) (1-p0)*dgamma(abs(x)+0.01,shape=sh,rate=ra)/fxinit(x),marginal=pxbetatag)
#                integral <- integral0 + integralnon0 
#                if(log(integral)>3) print(c(grrow,tag,i,log(integral)))
#                return(log(integral))
#                })
#            }
            if(modus=="mixt") {
                resall <- sapply(1:ngrid,function(grrow){
                #grrow <-180
                pminus <- gridd[grrow,1]
                p0 <- gridd[grrow,2]
                mu <- gridd[grrow,3]
                stdev <- gridd[grrow,4]
                if(!zerofit) integral0 <- (pxbeta_eq0tag/f0init)*(p0) else integral0 <- exp(mlik0minmlik1[tag])*p0
                integralnon0 <- myinla.expectation(function(x) (1-p0)*(pminus*dnorm(x,mean=-mu,sd=stdev) + 
                (1-pminus)*dnorm(x,mean=mu,sd=stdev))/fxinit(x),marginal=pxbetatag)  
                integral <- integral0 + integralnon0 
                if(log(integral)>3) print(c(grrow,tag,i,log(integral)))
                return(log(integral))
                })
            }     
            logliksum <- logliksum+resall
            #loglikall <- c(loglikall,c(tagsinv[tag],(which(arr[tag,]==1))[i],resall))
        }
        return(logliksum)
    }
set.seed(34523)
wh <- sample(1:ntag,ntagcur)
#if(ncpus==1) 
reslogliks <- sapply(wh,tagfun) 
#else { mysfExport(forceexport=c("pxbeta","pxbeta_eq0"),forceexclude=c("fitall","fitall0"))
#print("Finished exporting")
#reslogliks <- sfSapply(wh,tagfun)}

proc.time()-pmt
total <- apply(reslogliks,1,sum)
ordt <- order(total,decreasing=TRUE)
print(cbind(gridd[ordt[1:10],],total[ordt[1:10]]))
if(ng < length(ntagreeks)) gridd <- gridd[ordt[1:ngridreeks[(ng+1)]],]
print(nrow(gridd))
}
   
total <- apply(reslogliks,1,sum)
ordt <- order(total,decreasing=TRUE)
allparam <- cbind(gridd[ordt,],sumloglik=total[ordt])
best <- allparam[1,]
if((refine == 1) & refinegrid){
print(paste("Current estimates",best))
if(modus!="mixt"){
    p0est <- best[1]
    vecsort <- sort(p0vec)
    wh <- which(vecsort==p0est)
    if(wh == length(vecsort)) vecmaxim <- 1 else vecmaxim <- vecsort[wh+1]
    if(wh == 1) vecminim <- 0.5 else vecminim <- vecsort[wh-1]
    p0vec <- seq(vecminim,vecmaxim,length.out=5)
    
    stest <- best[2]
    vecsort <- sort(sdvec)
    wh <- which(vecsort==stest)
    if(wh == length(vecsort)) vecmaxim <- 2*stest else vecmaxim <- vecsort[wh+1]
    if(wh == 1) vecminim <- 0.5*stest else vecminim <- vecsort[wh-1]
    sdvec <- seq(vecminim,vecmaxim,length.out=5)  
} else {
    pminest <- best[1]
    p0est <- best[2]
    muest <- best[3]
    stest <- best[4]
    
    vecsort <- sort(pminvec)
    wh <- which(vecsort==pminest)
    if(wh == length(vecsort)) vecmaxim <- 1 else vecmaxim <- vecsort[wh+1]
    if(wh == 1) vecminim <- 0 else vecminim <- vecsort[wh-1]
    if(!symmetric) pminvec <- seq(vecminim,vecmaxim,length.out=7) else pminvec <- 0.5

    vecsort <- sort(p0vec)
    wh <- which(vecsort==p0est)
    if(wh == length(vecsort)) vecmaxim <- 1 else vecmaxim <- vecsort[wh+1]
    if(wh == 1) vecminim <- 0.5 else vecminim <- vecsort[wh-1]
    p0vec <- seq(vecminim,vecmaxim,length.out=7)

    vecsort <- sort(sdvec)
    wh <- which(vecsort==stest)
    if(wh == length(vecsort)) vecmaxim <- 1.5*stest else vecmaxim <- vecsort[wh+1]
    if(wh == 1) vecminim <- 0.5*stest else vecminim <- vecsort[wh-1]
    sdvec <- seq(vecminim,vecmaxim,length.out=4)
    
    vecsort <- sort(meanvec)
    wh <- which(vecsort==muest)
    if(wh == length(vecsort)) vecmaxim <- 1.5*muest else vecmaxim <- vecsort[wh+1]
    if(wh == 1) vecminim <- 0.5*muest else vecminim <- vecsort[wh-1]
    meanvec <- seq(vecminim,vecmaxim,length.out=4)
}
print("Start refinement")
}
}
print(paste("Final estimates",best))
return(list(allparam=allparam,inputpar=inputpar,best=best))
} #END FUNCTION
