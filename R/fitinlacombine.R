fitinlacombine <-
function(fitlist,modus, probs=c(1,0),para = NULL, nsam, curvepred=NULL, safemode=TRUE){ #provides a sample from the empirical mixture of posteriors of the parameter to be shrunken; 
#fitlist: list of FitInlaAll objects from different mixture components
#modus: "fixed", "fixedasrandom", "random", "p0" or "disp" only one mixture allowed
#para: parameter name for which to combine the posteriors, only needed for fixed or random regression parameters
#curvepred: prediction from loess curve for log-dispersion
#safemode: if TRUE, only use tags for which both model provide a fit, in case of a mixture
#fitlist <- fit1;modus<- "random"; para = "timefac";nsam=10
#fitlist <- list(fitzinb2[[1]],fitzinb1[[1]]);probs=c(0.2,0.8);para="tissue";nsam=100;curvepred=NULL;safemode=TRUE;modus<- "random"
#fitlist = list(fitall);probs=diracprob; modus="disp";nsam=nsamtagfixed;safemode=safemode;nsam=100
if(modus=="fixed") {if(!is.null(para)){if(is.factor(try(get(para),silent=TRUE))) para <- fact2vec(para)}}

lf <- length(fitlist)
postparall <- c()
lpwall <- c()
for (i in 1:lf){
#i<-1
    fi <- fitlist[[i]]
    if(modus=="p0") {
        postparcomp <- unlist(lapply(fi, function(ex) 
                            {
                            #ex<-fi[[31]]
                            mhs <- ex$internal.marginals.hyperpar$"intern zero-probability" #logit p
                            if(is.null(mhs)) toret <- rep(-Inf,nsam) else {
                                #wh0 <- which(substr(rownames(ex$internal.summary.hyperpar),1,3)=="int")
                                wh0 <- which(substr(names(ex$internal.marginals.hyperpar),1,3)=="int")  
                                mss <- ex$internal.summary.hyper$mean[wh0]  
                                if(is.null(mss)) mss <- ex$internal.summary.hyper[wh0,1]
                                if(mss <= -3) toret <- rep(-Inf,nsam) else toret <- samplepost(mhs,nsam)
                                }
                            return(toret)
                            })) 
        para <- "p0" #just to make sure length(para)==1        
    }  
    if(modus=="disp") { 
        postparcomp <- lapply(fi, function(ex2) {
                            #ex2<-fi[[31]]
                            mhs <- ex2$internal.marginals.hyperpar$"log size" #log-size
                            if(is.null(mhs)) toret <- rep(Inf,nsam) else {
                                wh0 <- which(substr(names(ex2$internal.marginals.hyperpar),1,8)=="log size")
                                mss <- ex2$internal.summary.hyper$mean[wh0]  
                                if(is.null(mss)) mss <- ex2$internal.summary.hyper[wh0,1]
                                if(mss > log(1000)) toret <- rep(Inf,nsam) else toret <- samplepost(mhs,nsam)
                                }
                            return(toret)
                            })
        if(!is.null(curvepred))  postparcomp <- {as.numeric(sapply(1:length(postparcomp),function(i) (postparcomp[[i]]-curvepred[i])))} else 
        {postparcomp <- unlist(postparcomp)}
        para <- "dispersion" #just to make sure length(para)==1       
    } 
    
    if(modus=="err") { 
    #fi<- fitgauss
        postparcomp <- lapply(fi, function(ex2) {
                            #ex2<-fi[[1]]
                            mhs <- ex2$internal.marginals.hyperpar$"Log precision for the Gaussian observations" #log-precision
                            if(is.null(mhs)) toret <- rep(Inf,nsam) else {
                                wh0 <- which(names(ex2$internal.marginals.hyperpar)=="Log precision for the Gaussian observations")
                                mss <- ex2$internal.summary.hyper$mean[wh0] 
                                if(is.null(mss)) mss <- ex2$internal.summary.hyper[wh0,1]
                                toret <- samplepost(mhs,nsam)
                                }
                            return(toret)
                            })
        postparcomp <- unlist(postparcomp)
        para <- "precerr" #just to make sure length(para)==1       
    } 
    
    if(modus=="fixed") { 

        postparcomp <- unlist(lapply(fi, function(ex2) {
                            mhs <- ex2$marginals.fixed
                            if(is.null(mhs)) toret <- rep(0,nsam*length(para)) else {
                                nms <- names(mhs) 
                                wh <- which(nms %in% as.vector(para))
                                if(length(wh)==0) toret <- rep(0,nsam*length(para)) else {
                                    mhswh <- mhs[wh]
                                    toret <- unlist(lapply(mhswh,samplepost,reps=nsam))
                                    }
                                }
                            return(toret)
                            }))
        
    postparcomp <- sapply(postparcomp, function(x) max(min(x,5),-5))    #x larger than 100 means variance smaller than 1/100.         
    } 
    
    if(modus=="random") { 
        precpara <- paste("Log precision for",para)
        nch <- nchar(precpara)
        postparcomp <- unlist(lapply(fi, function(ex2) {
                            #ex2 <- fi[[1]]
                            mhs <- ex2$internal.marginals.hyperpar
                            if(is.null(mhs)) toret <- rep(Inf,nsam) else {
                                nms <- names(mhs) 
                                wh <- which(nms==precpara)
                                whl <- length(wh)
                                if(whl==0) {
                                nms2 <- sapply(nms,substr,1,nch)
                                wh <- which(nms2==precpara)
                                }
                                if(length(wh)==0) toret <- rep(Inf,nsam) else {
                                mhswh <- mhs[[wh]]
                                toret <- samplepost(mhswh,nsam)
                                }
                            }
                            return(toret)
                            }))    
        postparcomp <- sapply(postparcomp, function(x) if(x<log(1000000)) max(x,log(1/1000000)) else Inf)    #x larger than 100 means variance smaller than 1/100. 
    }  
    marglik <- unlist(lapply(fi,function(ex) {ml <- ex$mlik[2]
    if(is.null(ml)) return(-Inf) else return(ml)
    }))
    lpw <- marglik + log(probs[i]+10^{-20})
    lpwall <- cbind(lpwall, lpw)
    postparall <- cbind(postparall,postparcomp)  
}  
    
    if(safemode) maxml <- apply(lpwall,1,min) else maxml <- apply(lpwall,1,max) #remove tags for which any (safemode) / all models have failed to fit
    whremove <- which(maxml==-Inf)
    print(paste("ntags not used: ",length(whremove)))
    if(length(whremove) >= (length(lpw)/2) & safemode) {
    maxml <- apply(lpwall,1,max)
    print("Many non-used tags: switching mode")
    whremove <- which(maxml==-Inf)
    print(paste("ntags not used: ",length(whremove)))
    }
    if(length(whremove)>0){
        lpwall <- lpwall[-whremove,] 
        whremovepost <- as.vector(sapply(whremove,function(x) ((x-1)*nsam*length(para) + 1):(x*nsam*length(para))))
        postparall <- postparall[-whremovepost,]
    }
    
if(lf==1) return(postparcomp[!is.na(postparcomp)]) else { #mixture prior present

    weightvec <- apply(lpwall,1,function(vec) {
                            mvec <- max(vec) #factorize the max out for stability reasons [exp(-large)]
                            logsum <- mvec + log(sum(exp(vec-mvec)))
                            return(rep(exp(vec-logsum),nsam*length(para)))
                            })
    weightvecmat <- matrix(weightvec,ncol=lf,byrow=T)
    
    sampost <- sapply(1:nrow(postparall),function(i){
                vec <- postparall[i,]
                prbs <- weightvecmat[i,] 
                as.numeric(sample(vec,size=1,prob=prbs))
                })
    return(sampost[!is.na(sampost)])
    }
}
