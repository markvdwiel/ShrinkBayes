CombinePosteriors <-
function(fullpost1,fullpost2,shrinksimul,para,includelincomb=TRUE,margcomb=c("marginals.fixed"),ncpus=2,ndigits=3,updateby=1000){ 
#ndigits: nr of digits used for storing the result; purely for memory efficiency reasons
#fullpost1=fitzip_shr;fullpost2=fitzinb_shr;shrinksimul=shrinksimul;para="groupfac";includelincomb=TRUE;margcomb=c("marginals.fixed");ncpus=2;ndigits=3;updateby=5000
priorp1=shrinksimul$pmlist$mixp[1]
priors <- fullpost1[[2]]
fullpost1 <- fullpost1[[1]]
fullpost2 <- fullpost2[[1]]

if(length(para)==1) {if(is.factor(get(para))) para <- fact2vec(para)}

mliks1 <- mliks(fullpost1) 
mliks2 <- mliks(fullpost2) 
whm <- which(!is.na(mliks1))[1]
ell <- length(fullpost1[[whm]]$marginals.lincomb.derived)
if(includelincomb & ell > 0) {
    margcomb <- c(margcomb,"marginals.lincomb.derived")
    layer2 <- list(para,NULL)
    } else {
    layer2 <- list(para)
    }
    
    
indexfun <- function(fp,ml){
    whmin <- which(!is.na(ml))[1]
    posti <- fp[[whmin]]
    nam <- names(posti)
    wh <- match(margcomb,nam)
    lenm <- length(margcomb)
    whlist <- lapply(1:lenm,function(x){
        whi<- wh[x]
        if(is.null(layer2[[x]])) toret <- 1:length(posti[[whi]]) else {
            lr <- layer2[[x]]
            nm <- names(posti[[whi]])
            toret <- match(lr,nm)
        }
        return(toret)
        })
    return(list(wh,whlist))
}
wh1 <- indexfun(fullpost1,mliks1)
wh2 <- indexfun(fullpost2,mliks2)
marglist <- function(fp,wh12) {
mlik <- fp$mlik[1,]
if(is.null(mlik)) return(NULL) else
unlist(lapply(1:length(wh12[[1]]),function(x) fp[[(wh12[[1]][x])]][wh12[[2]][[x]]]),recursive=FALSE)
}

 
el <- length(fullpost1)

fullpost1 <- lapply(fullpost1,marglist,wh12=wh1) #overwrite to save memory space
fullpost2 <- lapply(fullpost2,marglist,wh12=wh2)

postmixprop <- function(i,ml1,ml2,pp1){
    ml1i <- ml1[i]
    ml2i<- ml2[i]
    if(is.na(ml1i)) p <- 0 
    if(is.na(ml2i)) p <- 1
    if(!is.na(ml1i) & !is.na(ml2i)) p <- 1/(1+((1-pp1)/pp1)*exp(ml2i-ml1i))
    return(p)
} 

allp1 <- sapply(1:length(mliks1),postmixprop,ml1=mliks1,ml2=mliks2,pp1=priorp1)
 
combineposteriors <- function(i, pl1, pl2, p1all,ndig=ndigits){ 
     #i<-7430;pl1<-fullpost1;pl2 <-fullpost2;p1all<- allp1
     #i is index, fullpost1 is list of inla objects using model 1, post idem model 2, p1: posterior probabilities model 1
     #if(is.wholenumber(i/500)) print(paste("Fitting tag nr ",i))
     fullpost1i <- pl1[[i]]
     fullpost2i <- pl2[[i]]
     p1 <- p1all[i]
     if(is.null(fullpost1i) & is.null(fullpost2i)) toret <-NULL
     if((is.null(fullpost1i) | p1<=0.01) & !is.null(fullpost2i)) toret <-fullpost2i
     if(!is.null(fullpost1i) & (is.null(fullpost2i) | p1>=0.99)) toret <- fullpost1i
     if(!is.null(fullpost1i) & p1>0.01 & p1<0.99 & !is.null(fullpost2i)) { 
     toret <- lapply(1:length(fullpost1i), function(j){ #j<-1
            marg10 <- try(myinla.smarginal(fullpost1i[[j]],len=150))
            marg20 <- try(myinla.smarginal(fullpost2i[[j]],len=150))
            if(class(marg10) != "try-error" & class(marg20) != "try-error"){
                marg1 <- cbind(marg10$x,marg10$y)
                marg2 <- cbind(marg20$x,marg20$y)
               
                xord <- order(c(marg1[,1],marg2[,1])) 
                 
                y1t <- myinla.dmarginal(marg2[,1],marg1,thrnrow=140)           
                y1 <- c(p1*marg1[,2],p1*y1t)
                y2t <- myinla.dmarginal(marg1[,1],marg2,thrnrow=140)
                y2 <- c((1-p1)*y2t,(1-p1)*marg2[,2])
                y <- y1+y2 
                
                x <- c(marg1[,1],marg2[,1])[xord]  
                y <- y[xord] 
                xy <- signif(cbind(x,y),ndig)
                if(sum(y1t)==0) xy <- marg2
                if(sum(y2t)==0) xy <- marg1
                return(xy)
                } else {
                if(class(marg10) == "try-error" & class(marg20) != "try-error") return(fullpost2i[[j]])
                if(class(marg20) == "try-error" & class(marg10) != "try-error") return(fullpost1i[[j]])
                if(class(marg10) == "try-error" & class(marg20) == "try-error") return(NULL)
                }
            })
        names(toret) <- names(fullpost1i)
        }
      return(toret)
}

if(ncpus==1) {res <- lapply(1:el,combineposteriors,pl1=fullpost1,pl2=fullpost2,p1all=allp1)} else {
    sfInit(parallel=TRUE,cpus=ncpus) 
    sfLibrary(INLA)
    sfLibrary(ShrinkBayes)
if(el<=updateby){
    print("Exporting to slaves")
    #sfExport("fullpost1","fullpost2","allp1","myinla.dmarginal","myinla.smarginal","linear.spline")
    sfExport("fullpost1","fullpost2","allp1")
    print("Exporting done, start combining posteriors")
    res <- sfLapply(1:el,combineposteriors,pl1=fullpost1,pl2=fullpost2,p1all=allp1)
    sfRemoveAll()
} else {
print(paste("Large objects...Job cut in pieces of",updateby))
        res <- as.list(rep(NA,el))
        sec <- seq(1,el,by=updateby)
        secl <- length(sec)
        if(sec[secl] > el-2) sec[secl] <- el+1 else sec <- c(sec,el+1)
        secl <- length(sec)
        for(k in 1:(secl-1)){
        fullpost1k <- fullpost1[(sec[k]:(sec[k+1]-1))]
        fullpost2k <- fullpost2[(sec[k]:(sec[k+1]-1))]
        allp1k <- allp1[(sec[k]:(sec[k+1]-1))]
        print("Exporting to slaves")
        #sfExport("fullpost1k","fullpost2k","allp1k","myinla.dmarginal","myinla.smarginal","linear.spline")
        sfExport("fullpost1k","fullpost2k","allp1k")
        print("Exporting done, start combining posteriors")
        resk <- sfLapply(1:(sec[k+1]-sec[k]),combineposteriors,pl1=fullpost1k,pl2=fullpost2k,p1all=allp1k)
        print(paste(sec[k+1]-1,"Cases done"))
        res[(sec[k]):(sec[k+1]-1)] <- resk
        sfRemoveAll()
        }
}
sfStop()
}
return(list(res=res,priors=priors))
}
