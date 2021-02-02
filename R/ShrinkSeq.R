ShrinkSeq <- function(form, dat, fams="zinb", 
shrinkfixed=NULL,shrinkaddfixed=NULL,shrinkrandom=NULL,shrinkaddrandom=NULL,shrinkdisp=TRUE,shrinkp0=FALSE,mixtdisp=FALSE,
mixtrand=FALSE,excludefornull=NULL,fixedmeanzero=FALSE, addfixedmeanzero=TRUE,curvedisp=FALSE, maxiter=15, ntag=ifelse(is.null(excludefornull),c(100,200,500,1000),c(1000)), ntagcurve=5000, fixed = c(0,1/10), 
addfixed = c(0,1/10), randomprec = c(1,10^(-5)),addrandomprec=c(1,10^(-5)), logdisp = c(0,0.01),diracprob0=ifelse((mixtrand | !is.null(excludefornull)), 0.8, 0.2),logitp0 = 
c(0,0.01), fixedseed = TRUE, ndraw = 10000, safemode=TRUE, tol=ifelse((mixtrand | !is.null(excludefornull)),0.005,0.01), 
tolrand = 0.02, mliktol = 0.1, designlist=NULL, ...) { 

##
#form = y ~ biopt2 + f(pers) + f(tissue) + deltapos + chemo;dat=mirnorm;shrinkfixed="biopt2";shrinkaddfixed=c("deltapos","chemo"); fams="zinb"; 
#shrinkrandom="tissue";shrinkdisp=TRUE;shrinkp0=FALSE;mixtdisp=FALSE;mixtrand=FALSE;excludefornull=c("deltapos","chemo");fixedmeanzero=FALSE; 
#addfixedmeanzero=TRUE;curvedisp=FALSE; maxiter=10;tol=0.01; tolrand = 0.02; ntag=c(100,200); ntagcurve=5000; fixed = c(0,1/10); 
#addfixed = c(0,1/10);randomprec = c(1,10^(-5)); logdisp = c(0,0.01);diracprob0=0.8;logitp0 = c(0,0.01); 
#fixedseed = TRUE; ndraw = 10000; safemode=TRUE;ncpus2use <- 8;spikeslab<-FALSE

#form = y ~ 1+g;shrinkfixed="g";shrinkaddfixed=NULL; dat=count;fams="zinb";
#shrinkrandom=NULL;shrinkdisp=TRUE;shrinkp0=FALSE;mixtdisp=TRUE;mixtrand=FALSE;fixedmeanzero=FALSE; 
#addfixedmeanzero=TRUE;curvedisp=FALSE; maxiter=10;tol=0.01; tolrand = 0.02; ntag=c(100,200); ntagcurve=5000; fixed = c(0,1/10); 
#addfixed = c(0,1/10);randomprec = c(1,10^(-5)); logdisp = c(0,0.01);diracprob0=0.8;logitp0 = c(0,0.01); 
#fixedseed = TRUE; ndraw = 10000; safemode=TRUE;ncpus2use <- 6;spikeslab<-FALSE;excludefornull=NULL

#form = y ~ 1 + PM + timepos + chemo + organ + f(indiv);shrinkfixed="PM";shrinkaddfixed=c("organ","chemo","timepos");
#dat=mirseqnorm;fams="zinb";shrinkaddrandom=NULL;addrandomprec <- c(1,10^(-5));
#shrinkrandom=NULL;shrinkdisp=TRUE;shrinkp0=FALSE;mixtdisp=FALSE;mixtrand=FALSE;fixedmeanzero=FALSE; 
#addfixedmeanzero=TRUE;curvedisp=FALSE; maxiter=10;tol=0.01; tolrand = 0.02; ntag=c(100,200); ntagcurve=5000; fixed = c(0,1/10); 
#addfixed = c(0,1/10);randomprec = c(1,10^(-5)); logdisp = c(0,0.01);diracprob0=0.8;logitp0 = c(0,0.01); 
#fixedseed = TRUE; ndraw = 10000; safemode=TRUE;ncpus2use <- 2;excludefornull="organ";mliktol=0.1


# 
# fams="zinb"; form = y ~ 1 + groupfac + batch + f(pers,model="iid"); dat=CAGEdata; shrinkaddfixed=NULL;shrinkaddrandom=NULL;shrinkdisp=TRUE;shrinkp0=FALSE;mixtdisp=FALSE;mixtrand=FALSE;excludefornull=NULL;fixedmeanzero=FALSE;
# addfixedmeanzero=TRUE;curvedisp=FALSE; maxiter=15; ntagcurve=5000; fixed = c(0,1/10); addfixed = c(0,1/10); randomprec =
# c(1,10^(-5)); addrandomprec=c(1,10^(-5));logdisp = c(0,0.01);diracprob0=0.8;logitp0 =
# c(0,0.01); fixedseed = TRUE; ndraw = 10000; safemode=TRUE; tol=0.005;
# tolrand = 0.02; mliktol = 0.1;
# maxiter=1;mixtrand=FALSE;ntag=c(50,100)
# #
# fams="zinb"; dat=cl.4; shrinkaddfixed=NULL;shrinkdisp=TRUE;shrinkp0=FALSE;mixtdisp=FALSE;mixtrand=FALSE;excludefornull=NULL;fixedmeanzero=FALSE;
# addfixedmeanzero=TRUE;curvedisp=TRUE; maxiter=15; ntagcurve=500; fixed = c(0,1/10); addfixed = c(0,1/10); randomprec =
# c(1,10^(-5)); logdisp = c(0,0.01);diracprob0=0.8;logitp0 =
# c(0,0.01); fixedseed = TRUE; ndraw = 10000; safemode=TRUE; tol=0.005;
# tolrand = 0.02; mliktol = 0.1;
  
#form=f.4; dat=cl.4;shrinkfixed="g.4";shrinkrandom=NULL;ncpus=8;
#maxiter=1;ntag=c(50,100)
#form=f.4, dat=cl.4, fams="nb", shrinkfixed="g.4", curvedisp=T, ncpus=4

#
#shrinksimul <- ShrinkSeq(form=form, dat=mirnorm,shrinkfixed="biopt2",shrinkrandom="tissue",shrinkaddfixed=c("chemo","deltapos"),
#fixedmeanzero=FALSE, addfixedmeanzero=FALSE, shrinkdisp=TRUE,mixtrand=TRUE,ncpus=ncpus2use,maxiter=2,ntags=c(100))


#form: see INLA manual. And examples therein 
#maxiter: maximum number of iterations per cycle
#tol: tolerance for KS distance. If KS distance for all priors is smaller than tol for two subsequent iterations, the algorithm stops
#tolrand: tolerance for KS distance for random effects parameter. See tol. Since posteriors of regression parameters are rather insensitive to the left-tail of this distribution we advise to use a larger value than tol.
#ntag: sequence of number of tags used to fit the priors. Very important for computing time
#shrinkfixed: which fixed parameter to shrink; this contains the main parameter of interest; or NULL. Can be a factor
#shrinkaddfixed: which fixed nuisance parameter to shrink. Can be a factor
#shrinkrandom: which random parameter to shrink; or NULL.Can be a factor.
#shrinkdisp: FALSE/TRUE: shrink overdispersion parameter?
#shrinkp0: shrink zero-inflation parameter?
#fixedmeanzero: if TRUE the Gaussian prior of the fixed effect is assumed to have mean 0
#mixtdisp: use mixture model (dirac plus Gaussian) for overdispersion prior?
#diracprob: mixture proportions for disp0 and disp positive
#curvedisp: should it apply nonparamteric curve fitting to shrink the logdisp?
#ntagcurve: number of tags used for initial curve fitting
#fixed: initial values for fixed shrinkage parameter; Gaussian, mu and prec
#randomprec: initial values for random shrinkage parameter; Gamma, alpha1 and alpha2
#logdisp: initial values for Gaussian (component of) prior for overdispersion; mu, prec
#diracprob: initial values for mixture proportions of dirac-Gauss mixture for overdispersion. For value is q0 (dirac)/
#logitp0: initial values for zero-inflation shrinkage parameter; Gaussian, mu and prec
#fixedseed: fix the seeds per component in ntag for reproducible results (same subsets of tags are used for given value of ntag)
#ndraw: number of draws from the posterior to be used to fit the hyperprior
#safemode: if TRUE, only use tags for which both model provide a fit, in case of a mixture

#form = y ~ 1 + grCN + f(batch1CN) + f(batch2CN) + f(batch3CN) + offset(offsetsCN);shrinkfixed="grCN";shrinkaddfixed=NULL;dat=dat;fams="zinb";
#shrinkrandom="batch1CN";shrinkaddrandom=c("batch2CN","batch3CN");shrinkdisp=TRUE;shrinkp0=FALSE;mixtdisp=FALSE;mixtrand=FALSE;fixedmeanzero=FALSE; 
#addfixedmeanzero=TRUE;curvedisp=FALSE; maxiter=2;tol=0.01; tolrand = 0.02; ntag=c(5,10); ntagcurve=5000; fixed = c(0,1/10); 
#addfixed = c(0,1/10);randomprec = c(1,10^(-5)); addrandomprec=c(1,10^(-5));logdisp = c(0,0.01);diracprob0=0.8;logitp0 = c(0,0.01); 
#fixedseed = TRUE; ndraw = 10000; safemode=TRUE;ncpus2use <- 1;excludefornull="grCN";mliktol=0.1


if(mixtdisp & mixtrand) {
    print("Only one mixture prior can be used: mixtdisp is set to FALSE")
    mixtdisp <- FALSE
    }
    
if(mixtdisp & !is.null(excludefornull)) {
    print("Only one mixture prior can be used: mixtdisp is set to FALSE")
    mixtdisp <- FALSE
    }
    
if(mixtrand & !is.null(excludefornull)) {
    print("Only one mixture prior can be used: excludefornull is set to NULL")
    exlcludefornull <- NULL
    }


if(!is.null(excludefornull)){
nexcl <- length(excludefornull)
form0 <- form
for(j in 1:nexcl) form0 <- replacerand0(form0,excludefornull[j])
}

#small auxiliary functions
repNA <- function(x) {if(is.na(x)) return(-10^10) else return(x)}
iselement2 <- function(elem,set) if(is.null(set) | is.null(elem)) return(FALSE) else return(is.element(elem,set))

#NEW SS > 1.8    
els <- length(shrinkaddfixed)
if(els >= 2){
if(mode(addfixedmeanzero) != "list") addfixedmeanzero <- as.list(rep(addfixedmeanzero,els))
if(mode(addfixed) != "list") {add0 <- addfixed; addfixed <- list(); for(i in 1:els) addfixed <- c(addfixed,list(add0))}
} else {
addfixedmeanzero <- list(addfixedmeanzero)
addfixed <- list(addfixed)
}

#NEW210    
elsr <- length(shrinkaddrandom)
if(elsr >= 2){
add0 <- addrandomprec; addrandom <- list(); for(i in 1:elsr) addrandom <- c(addrandom,list(add0))
} else {
addrandom <- list(addrandomprec)
}

formch<-deparse(form)

inputpar <- list(formch,shrinkfixed, shrinkaddfixed,shrinkrandom,shrinkaddrandom,shrinkdisp,shrinkp0,mixtdisp,mixtrand,excludefornull,fixedmeanzero, 
addfixedmeanzero,curvedisp, maxiter,tol, tolrand, ntagcurve, fixed, addfixed, randomprec, addrandomprec, 
logdisp,diracprob0,logitp0, fixedseed, ndraw, safemode, mliktol) 
names(inputpar) <- c("form", "shrinkfixed", 
"shrinkaddfixed", "shrinkrandom","shrinkaddrandom", "shrinkdisp", "shrinkp0", "mixtdisp", "mixtrand","excludefornull", "fixedmeanzero", 
"addfixedmeanzero", "curvedisp", "maxiter", "tol", "tolrand", "ntagcurve", "fixed", "addfixed", "randomprec", "addrandomprec",
"logdisp", "diracprob0", "logitp0", "fixedseed", "ndraw", "safemode", "mliktol") 

diracprob <- c(diracprob0,1-diracprob0)


#### here starts function
if(class(dat)[1] =="list")  ntagtot <- ntagtotal <- length(dat) else ntagtot <- ntagtotal <- nrow(dat) #NEW


if (length(fams)==1) fams <- rep(fams,ntagtot)


#ntagtotal <- nrow(dat)
wh2large <- which(ntag >= ntagtotal)
if(length(wh2large)>0){
ntag <- c(ntag[-wh2large],ntagtotal)
}
lngene<-length(ntag)

#### CREATE INITIAL PARAMETER LIST  #####
pmlist <- list()
#if(!is.null(shrinkfixed)) {

        pmlist <- c(pmlist,list(mufixed=fixed[1],precfixed=fixed[2]))
#}

if(!is.null(shrinkaddfixed)) {
#NEW >= 1.8
        for(i in 1:els) {addf <- as.numeric(addfixed[[i]]);name <- c(paste("mu",shrinkaddfixed[i],sep=""),
        paste("prec",shrinkaddfixed[i],sep="")); toadd <- list(muaddfixed=addf[1],precaddfixed=addf[2]);
        names(toadd) <- name;
        pmlist <- c(pmlist,toadd)}
} else {
pmlist <- c(pmlist,list(muaddfixed=addfixed[[1]][1],precaddfixed=addfixed[[1]][2]))
}



#if(!is.null(shrinkrandom){
        pmlist <- c(pmlist,list(shaperand = randomprec[1],raterand=randomprec[2]))
#}

#NEW210
elsr <- length(shrinkaddrandom)
if(!is.null(shrinkaddrandom)) {
        for(i in 1:elsr) {addr <- addrandomprec;name <- c(paste("shape",shrinkaddrandom[i],sep=""),
        paste("rate",shrinkaddrandom[i],sep="")); toadd <- list(shapeaddr=addr[1],rateaddr=addr[2]);
        names(toadd) <- name;
        pmlist <- c(pmlist,toadd)}
} else {
pmlist <- c(pmlist,list(shapeaddr=addrandomprec[1],rateaddr=addrandomprec[2]))
}
    
#if(shrinkp0){
        pmlist <- c(pmlist,list(mup0=logitp0[1],precp0=logitp0[2]))
 #}   
# if(shrinkdisp){
        pmlist <- c(pmlist,list(mudisp=logdisp[1],precdisp=logdisp[2],mixp=diracprob))

print(pmlist)   
ntagused <- NA 
paraprev <- c(unlist(pmlist),nfeat=NA, meanmlik = NA)
paraall <- paraprev 
ksall <- c()

cf <- list(mean=0, prec = 0.01) #default prior for fixed effect, NO SHRINKAGE
          

 ###### INITIAL FIT WHEN CURVE FITTING FOR LOGDISP IS REQUIRED  #########
if(curvedisp){
print(paste("Fitting curve based on",ntagcurve,"observations"))
    ngenej <- ntagcurve
    famlevels <- unique(fams)
    fllen <- length(famlevels)
    if(fixedseed){set.seed(43545126)}  #prefer to use a fixed seed per loop to have reproducible results
    if(fllen==1){ #when fams differ make sure all fams are correctly represented
        sel <- sample(1:ntagtotal,min(ntagtotal,ngenej))
        } else {
            sel <- c()
            for(j1 in 1:fllen){
            wh <- which(fams==famlevels[j1])
            prop <- length(wh)/ntagtotal
            nsample <- min(length(wh),round(prop*ngenej))
            sel <- c(sel,sample(wh,nsample))
            }
        }
    if(class(dat)[1] =="list") datshrinkj <- dat[sel] else datshrinkj <- dat[sel,]   #NEW
    designlistj <- designlist[sel] #NEW
    famsj <- fams[sel]
    #fitallinit <- FitInlaAll(form,datshrinkj,famsj, logdisp, logitp0, cf=cf, control.compute=list(dic=F, mlik=T, cpo=F),curvedispfun=NULL)
    fitallinit <- FitInlaAll(form,datshrinkj,famsj, logdisp=logdisp, logitp0=logitp0, cf=cf, control.compute=list(dic=F, mlik=T, cpo=F),curvedispfun=NULL, designlist=designlistj,...)
    logsums <- log(apply(datshrinkj,1,sum))
    logdisps <- unlist(lapply(fitallinit,function(ex) {
    #ex <- fitallinit[[1]]
            if(!is.null(ex$internal.summary.hyperpar)){
            wh0 <- which(substr(rownames(ex$internal.summary.hyperpar),1,8)=="log size")
            mss <- ex$internal.summary.hyper[wh0,1]
            } else mss <- NA
            return(mss)
            }))
    whichna <- which(is.na(logdisps))
    if(length(whichna)>0){
    logdisps <- logdisps[-whichna] 
    logsums <- logsums[-whichna]  
    }                
    predLoess <- function(logDisp, logTotalCounts){
    fit <- loess(logDisp~logTotalCounts, span = 0.75, degree=2)
    fun = function(x) {
        if(x<= min(logsums)) x <- min(logsums) 
        if(x>= max(logsums)) x <- max(logsums)
        predict(fit,x)
    }
    return(fun)
    }
    curvedispfun <- predLoess(logdisps,logsums)
} else {
    curvedispfun <- NULL
} 
     

###### START FOR-LOOP OVER NTAG INCLUDED (OUTER-LOOP_  ##########################

for(j in 1:lngene){
#j<-1
    ngenej <- ntag[j]
    iter <- 1
    change<-1
    moreiter <- TRUE
    famlevels <- unique(fams)
    fllen <- length(famlevels)
    if(fixedseed){set.seed(43545126+j)}  #prefer to use a fixed seed per loop to have reproducible results
    if(fllen==1){ #when fams differ make sure all fams are correctly represented
        sel <- sample(1:ntagtotal,min(ntagtotal,ngenej))
        } else {
            sel <- c()
            for(j1 in 1:fllen){
            wh <- which(fams==famlevels[j1])
            prop <- length(wh)/ntagtotal
            nsample <- min(length(wh),round(prop*ngenej))
            sel <- c(sel,sample(wh,nsample))
            }
        }
    if(class(dat)[1] =="list") datshrinkj <- dat[sel] else datshrinkj <- dat[sel,]   #NEW
    designlistj <- designlist[sel] #NEW
    famsj <- fams[sel]
        
    famsjdisp0 <- sapply(famsj,function(fam) if(fam=="zinb" | fam=="zip") "zip" else "poisson") #families when disp pm set to 0
    famsjp00 <- sapply(famsj,function(fam) if(fam=="zinb" | fam=="nb") "nb" else "poisson") ##families when p0 pm set to 0
    
    ###### START INNER-LOOP UNTIL CONVERGENCE OR MAXITER  ##########################
    mlikmeanall <- c() 
    mlikprev <- -Inf
    while(moreiter){
        print(paste("iter=",iter))
        print(paste("ntagsused=",ntag[j]))
        if(diracprob[1]==0) {
            mixtrand<-FALSE;mixtdisp<- FALSE; excludefornull <- NULL; 
            print("Point mass has zero mass...fitting full model only.")
        }
            ########### FIRST UPDATE INLA CALL #################
        #NEW >= 1.8 
            cf <- list(mean=list(default=0), prec = list(default=0.01))
            if(!is.null(shrinkfixed)){
                if(is.factor(try(get(shrinkfixed),silent=TRUE))) shrinkfixed <- fact2vec(shrinkfixed)
                   nm <- names(cf[[1]])
                   for(j1 in 1:length(shrinkfixed)){       
                    cf$mean <- c(cf$mean,list(fixed[1]))
                    cf$prec <- c(cf$prec,list(fixed[2]))
                    }
                    names(cf[[1]]) <- c(nm,shrinkfixed)
                    names(cf[[2]]) <- c(nm,shrinkfixed) 
                    } 

            if(!is.null(shrinkaddfixed)){ 
                for(i in 1:els){
                #i<-1
                    shrinkaddfixedi <- shrinkaddfixed[i]
                    addfixedi <- as.numeric(addfixed[[i]])   
                    nm <- names(cf[[1]])       
                    if(is.factor(try(get(shrinkaddfixedi),silent=TRUE))) shrinkaddfixedi <- fact2vec(shrinkaddfixedi)
                    for(j1 in 1:length(shrinkaddfixedi)){       
                    cf$mean <- c(cf$mean,list(addfixedi[1]))
                    cf$prec <- c(cf$prec,list(addfixedi[2]))
                    }
                    names(cf[[1]]) <- c(nm,shrinkaddfixedi)
                    names(cf[[2]]) <- c(nm,shrinkaddfixedi)  
                    }
                }

            
            if(!is.null(shrinkrandom)){
                form <- randreplace(form,shrinkrandom,randomprec)
                formrand0 <- replacerand0(form,shrinkrandom)
                if(!is.null(excludefornull) & !iselement2(inputpar$shrinkrandom,excludefornull)){
                form0 <- randreplace(form0,shrinkrandom,randomprec)
                }
                }
                
            #NEW210
            if(!is.null(shrinkaddrandom)){
                for(i in 1:elsr){
                form <- randreplace(form,shrinkaddrandom[i],addrandom[[i]])
                if(!is.null(excludefornull) & !iselement2(inputpar$shrinkaddrandom[i],excludefornull)){
                form0 <- randreplace(form0,shrinkaddrandom[i],addrandom[[i]])
                }
                }
                }
            ################# END UPDATE INLA CALL ############################################
            
            ################# APPLY INLA FOR THE CURRENT SUBSET OF TAGS #######################
        
        #### !!!!!!!!!!! INSERT DOT NOTATION BELOW !!!!!!!!!!  #####
            
           
            nsamtag <- ceiling(ndraw/ngenej) #number of samples per posterior
            nsamtagfixed <- ceiling(ndraw/(max(1,length(shrinkfixed))*ngenej))
            

            
           fitall <- FitInlaAll(form,datshrinkj,famsj, logdisp=logdisp, curvedispfun=curvedispfun, logitp0=logitp0, cf=cf, control.compute=list(dic=F, mlik=T, cpo=F), ndigits=5, designlist=designlistj, ...)
           #fitall <- FitInlaAll(form,datshrinkj,famsj, logdisp=logdisp, curvedispfun=curvedispfun, logitp0=logitp0, control.compute=list(dic=F, mlik=T, cpo=F),cf=cf, ncpus=ncpus2use, ndigits=5)
           mliksall <- mliks(fitall)
            if(mixtdisp) {
                fitalldisp0 <- FitInlaAll(form,datshrinkj,famsjdisp0, logdisp=logdisp,curvedispfun=curvedispfun, logitp0=logitp0, control.compute=list(dic=F, mlik=T, cpo=F),cf=cf,ndigits=5, designlist=designlistj,...)         
                #fitalldisp0 <- FitInlaAll(form,datshrinkj,famsjdisp0, logdisp, curvedispfun=curvedispfun, logitp0, cf=cf, control.compute=list(dic=F, mlik=T, cpo=F), ndigits=5, ncpus=ncpus2use)
                mliksall0 <- mliks(fitalldisp0)       
                if(!is.null(shrinkfixed)) postfixed <- fitinlacombine(list(fitalldisp0,fitall),probs=diracprob, modus="fixed",para=shrinkfixed,nsam=nsamtagfixed,safemode=safemode)
                if(!is.null(shrinkaddfixed)){
                    postaddfixed <- list()
                    for(i in 1:els){
                    shrinkaddfixedi <- shrinkaddfixed[i]
                    postaddfixedi <- fitinlacombine(list(fitalldisp0,fitall),probs=diracprob, modus="fixed",para=shrinkaddfixedi,nsam=nsamtagfixed,safemode=safemode)
                    postaddfixed <- c(postaddfixed,list(postaddfixedi))
                    }
                }
                if(!is.null(shrinkrandom)) postrandom <- fitinlacombine(list(fitalldisp0,fitall),probs=diracprob, modus="random",para=shrinkrandom,nsam=nsamtag,safemode=safemode)
                
                #NEW210
                if(!is.null(shrinkaddrandom)){
                    postaddrandom <- list()
                    for(i in 1:elsr){
                    shrinkaddrandomi <- shrinkaddrandom[i]
                    postaddrandomi <- fitinlacombine(list(fitalldisp0,fitall),probs=diracprob, modus="random",para=shrinkaddrandomi,nsam=nsamtag,safemode=safemode)
                    postaddrandom <- c(postaddrandom,list(postaddrandomi))
                    }
                }
                
                
                if(shrinkp0) postp0 <- fitinlacombine(list(fitalldisp0,fitall),probs=diracprob, modus="p0",nsam=nsamtag,safemode=safemode) #PROBS IS NOT RELEVANT IN THIS CASE
                if(shrinkdisp) postdisp <- fitinlacombine(list(fitalldisp0,fitall),probs=diracprob, modus="disp",nsam=nsamtag,safemode=safemode)
            } #END mixtdisp
            
            if(mixtrand) {
                fitallrand0 <- FitInlaAll(formrand0,datshrinkj, logdisp=logdisp,curvedispfun=curvedispfun, logitp0=logitp0, control.compute=list(dic=F, mlik=T, cpo=F),cf=cf, ndigits=5,designlist=designlistj, ...)      
                mliksall0 <- mliks(fitallrand0) 
                #fitallrand0 <- FitInlaAll(formrand0,datshrinkj, logdisp=logdisp,curvedispfun=curvedispfun, logitp0=logitp0, control.compute=list(dic=F, mlik=T, cpo=F),cf=cf, ndigits=5)  
                if(!is.null(shrinkfixed)) postfixed <- fitinlacombine(list(fitallrand0,fitall),probs=diracprob, modus="fixed",para=shrinkfixed,nsam=nsamtagfixed,safemode=safemode)
                if(!is.null(shrinkaddfixed)) {
                    postaddfixed <- list()
                    for(i in 1:els){
                    shrinkaddfixedi <- shrinkaddfixed[[i]]
                    postaddfixedi <- fitinlacombine(list(fitallrand0,fitall),probs=diracprob, modus="fixed",para=shrinkaddfixedi,nsam=nsamtagfixed,safemode=safemode)
                    postaddfixed <- c(postaddfixed,list(postaddfixedi))
                    }
                }  
                if(!is.null(shrinkrandom)) postrandom <- fitinlacombine(list(fitallrand0,fitall),probs=diracprob, modus="random",para=shrinkrandom,nsam=nsamtag,safemode=safemode)
                
                #NEW210
                if(!is.null(shrinkaddrandom)){
                    postaddrandom <- list()
                    for(i in 1:elsr){
                    shrinkaddrandomi <- shrinkaddrandom[i]
                    postaddrandomi <- fitinlacombine(list(fitallrand0,fitall),probs=diracprob, modus="random",para=shrinkaddrandomi,nsam=nsamtag,safemode=safemode)
                    postaddrandom <- c(postaddrandom,list(postaddrandomi))
                    }
                }
                
                if(shrinkp0) postp0 <- fitinlacombine(list(fitallrand0,fitall),probs=diracprob, modus="p0",nsam=nsamtag,safemode=safemode) #PROBS IS NOT RELEVANT IN THIS CASE
                if(shrinkdisp) if(is.null(curvedispfun)) {
                postdisp <- fitinlacombine(list(fitallrand0,fitall),probs=diracprob, modus="disp",nsam=nsamtag,safemode=safemode)
                } else {
                curvepred <- sapply(log(apply(datshrinkj,1,sum)),curvedispfun)
                postdisp <- fitinlacombine(list(fitallrand0,fitall),probs=diracprob, modus="disp",nsam=nsamtag,curvepred=curvepred,safemode=safemode)
                }        
            }
            
            if(!is.null(excludefornull)) {
                fitall0 <- FitInlaAll(form0,datshrinkj, logdisp=logdisp,curvedispfun=curvedispfun, logitp0=logitp0, control.compute=list(dic=F, mlik=T, cpo=F),cf=cf,ndigits=5,designlist=designlistj,...)      
                #fitall0 <- FitInlaAll(form0,datshrinkj, logdisp=logdisp,curvedispfun=curvedispfun, logitp0=logitp0, control.compute=list(dic=F, mlik=T, cpo=F),cf=cf,ndigits=5)  
                mliksall0 <- mliks(fitall0)
                if(!is.null(shrinkfixed)) postfixed <- fitinlacombine(list(fitall0,fitall),probs=diracprob, modus="fixed",para=shrinkfixed,nsam=nsamtagfixed,safemode=safemode)
                if(!is.null(shrinkaddfixed)) {
                    postaddfixed <- list()
                    for(i in 1:els){
                    shrinkaddfixedi <- shrinkaddfixed[[i]]
                    postaddfixedi <- fitinlacombine(list(fitall0,fitall),probs=diracprob, modus="fixed",para=shrinkaddfixedi,nsam=nsamtagfixed,safemode=safemode)
                    postaddfixed <- c(postaddfixed,list(postaddfixedi))
                    }
                }  
                if(!is.null(shrinkrandom)) postrandom <- fitinlacombine(list(fitall0,fitall),probs=diracprob, modus="random",para=shrinkrandom,nsam=nsamtag,safemode=safemode)
                
                #NEW210
                if(!is.null(shrinkaddrandom)){
                    postaddrandom <- list()
                    for(i in 1:elsr){
                    shrinkaddrandomi <- shrinkaddrandom[i]
                    postaddrandomi <- fitinlacombine(list(fitall0,fitall),probs=diracprob, modus="random",para=shrinkaddrandomi,nsam=nsamtag,safemode=safemode)
                    postaddrandom <- c(postaddrandom,list(postaddrandomi))
                    }
                }
                
                if(shrinkp0) postp0 <- fitinlacombine(list(fitall0,fitall),probs=diracprob, modus="p0",nsam=nsamtag,safemode=safemode) #PROBS IS NOT RELEVANT IN THIS CASE
                if(shrinkdisp) if(is.null(curvedispfun)) {
                postdisp <- fitinlacombine(list(fitall0,fitall),probs=diracprob, modus="disp",nsam=nsamtag,safemode=safemode)
                } else {
                curvepred <- sapply(log(apply(datshrinkj,1,sum)),curvedispfun)
                postdisp <- fitinlacombine(list(fitall0,fitall),probs=diracprob, modus="disp",nsam=nsamtag,curvepred=curvepred,safemode=safemode)
                }        
            } 
            
            if(!mixtrand & !mixtdisp & is.null(excludefornull)){
                if(!is.null(shrinkfixed)) postfixed <- fitinlacombine(list(fitall), modus="fixed",para=shrinkfixed,nsam=nsamtagfixed)
                if(!is.null(shrinkaddfixed)) {
                        postaddfixed <- list()
                        for(i in 1:els){
                        shrinkaddfixedi <- shrinkaddfixed[i]
                        postaddfixedi <- fitinlacombine(list(fitall), modus="fixed",para=shrinkaddfixedi,nsam=nsamtagfixed,safemode=safemode)
                        postaddfixed <- c(postaddfixed,list(postaddfixedi))
                        }
                    }                  
                if(!is.null(shrinkrandom)) postrandom <- fitinlacombine(list(fitall), modus="random",para=shrinkrandom,nsam=nsamtag)
                
                #NEW210
                if(!is.null(shrinkaddrandom)){
                    postaddrandom <- list()
                    for(i in 1:elsr){
                    shrinkaddrandomi <- shrinkaddrandom[i]
                    postaddrandomi <- fitinlacombine(list(fitall),probs=diracprob, modus="random",para=shrinkaddrandomi,nsam=nsamtag,safemode=safemode)
                    postaddrandom <- c(postaddrandom,list(postaddrandomi))
                    }
                }
                if(shrinkp0) postp0 <- fitinlacombine(list(fitall), modus="p0",nsam=nsamtag) #PROBS IS NOT RELEVANT IN THIS CASE
                if(shrinkdisp) if(is.null(curvedispfun)) {
                    postdisp <- fitinlacombine(list(fitall), modus="disp",nsam=nsamtag)
                    } else {
                    curvepred <- sapply(log(apply(datshrinkj,1,sum)),curvedispfun)
                    postdisp <- fitinlacombine(list(fitall), modus="disp",nsam=nsamtag,curvepred=curvepred)
                    }
                }
           
           if(!mixtrand & !mixtdisp & is.null(excludefornull)) mlikmean <- mean(mliksall,na.rm=T) else {
                whna <- union(which(is.na(mliksall)),which(is.na(mliksall0)))
                if(length(whna)>0) {mliksall <- mliksall[-whna]; mliksall0 <- mliksall0[-whna]}
                #mlikmean <- mean(log(diracprob[1]) + mliksall0 + log(1+exp((log(diracprob[2])+mliksall-log(diracprob[1])-mliksall0))),na.rm=T)
                if(diracprob[1]==0) mlikmean <- mean(mliksall,na.rm=T)
                if(diracprob[2]==0) mlikmean <- mean(mliksall0,na.rm=T)
                if(diracprob[1]!=0 & diracprob[2]!=0) mlikmean <- mean(log(diracprob[2]) + mliksall + log(1+exp((log(diracprob[1])+mliksall0-log(diracprob[2])-mliksall))),na.rm=T)
                }
           mlikmeanall <- c(mlikmeanall,mlikmean)
           print(mlikmeanall)
           # plot(density(postfixed))
            ############### END APPLY INLA ######################################################
          
            ################### ESTIMATE HYPERPRIOR PARAMETERS FOR NON-MIXTURE PRIORS ##################################
          mlikrel <- (mlikmean-mlikprev)/abs(mlikmean)
          print(mlikrel)
          mlikconv <-  (mlikrel <= (mliktol/100))
          
          #NEW210  
          if(mlikconv){
            fixed <- fixedprev
            addfixed <- addfixedprev
            addrandom <- addrandomprev
            randomprec <- randomprecprev
            logitp0<- logitp0prev 
            logdisp <- logdispprev
            diracprob <- diracprobprev
            } else {  
            fixedprev <- fixed
            addfixedprev <- addfixed
            addrandomprev <- addrandom
            randomprecprev <- randomprec
            logitp0prev <- logitp0
            logdispprev <- logdisp
            diracprobprev <- diracprob
            
            #NEW210
            if(is.null(excludefornull)) p0 <- 0 else {
            mlik <- mliks(fitall)
            mlik0 <- mliks(fitall0)
            repNA <- function(x) {if(is.na(x)) return(-10^10) else return(x)}
            mlik <- sapply(mlik,repNA)
            mlik0 <- sapply(mlik0,repNA)
            
            #added for numerical stability
            mlik0 <- mlik + sapply(mlik0 - mlik,function(x) min(40,max(x,-40)))
            
            maxlik <- as.numeric(apply(cbind(mlik,mlik0),1,max))
            p0start <- length(which((mlik0-mlik)>0))/length(mlik)
            
            liktot <- function(p0=0.5) {-sum(log(p0*exp(mlik0-maxlik) + (1-p0)*exp(mlik-maxlik)))}
            #res1 <- mle(liktot,start=list(p0=p0start),lower=c(0.001),upper=c(0.999),hessian=FALSE,method="L-BFGS-B")
            res2 <- optimize(liktot,lower=0,upper=1, maximum = FALSE) #seems more stable than mle
            p0 <- res2$minimum        
            }
            
               
             if(!is.null(shrinkfixed)) {
                if(is.element(inputpar$shrinkfixed,excludefornull)) {
                    wh0 <- which(postfixed==0)
                    p0 <- length(wh0)/length(postfixed)
                    postfixed <- postfixed[-wh0]
                    }
                if(length(postfixed)>1){   
                if(fixedmeanzero) fixed <- gaussML(postfixed) else fixed <- gaussMLwithmean(postfixed)  #gaussian prior
                }
            }
            
             if(!is.null(shrinkaddfixed)) {
        #NEW SS >= 1.8
                addfixed <- list()
                for(i in 1:els){
                #i<-1
                addfixedmeanzeroi <- addfixedmeanzero[[i]]
                postaddfixedi <- unlist(postaddfixed[[i]])
                if(is.element(inputpar$shrinkaddfixed[i],excludefornull)) {
                    wh0 <- which(postaddfixedi==0)
                    p0 <- length(wh0)/length(postaddfixedi)
                    postaddfixedi <- postaddfixedi[-wh0]
                    }
                if(length(postaddfixedi)>1){   
                if(addfixedmeanzeroi) addfixedi <- gaussML(postaddfixedi) else addfixedi <- gaussMLwithmean(postaddfixedi)
                }
                addfixed <- c(addfixed,list(addfixedi))
                }
            }
            
            if(!is.null(shrinkrandom)){
            if(mixtrand){
                   mlm <- mixture2loggammaML(postrandom)
                   randomprec <- c(mlm[2],mlm[3])
                   p0 <- mlm[1]
            } else {
                 if(iselement2(inputpar$shrinkrandom,excludefornull)) {
                 wh0 <- which(postrandom==0)
                 p0 <- length(wh0)/length(postrandom)
                 postrandom <- postrandom[-wh0]
                 }
                 if(length(postrandom)>1){   
                  randomprec <- loggammaML(postrandom)
                  }
                  }
            }
        #NEW210
             if(!is.null(shrinkaddrandom)) {
                addrandom <- list()
                for(i in 1:elsr){
                 postaddrandomi <- unlist(postaddrandom[[i]])
                if(is.element(inputpar$shrinkaddrandom[i],excludefornull)) {
                    wh0 <- which(postaddrandomi==0)
                    p0 <- length(wh0)/length(postaddrandomi)
                    postaddrandomi <- postaddrandomi[-wh0]
                    }
                if(length(postaddrandomi)>1){   
                addrandomi <- loggammaML(postaddrandomi)
                } else addrandomi <- addrandomprec
                addrandom <- c(addrandom,list(addrandomi))
                }
            }
            
                
            if(shrinkp0){ 
                     logitp0 <- gaussMLwithmean(postp0) 
             }  
          
             
             
            if(shrinkdisp){
                if(mixtdisp){  
                    mixfit <- mixture2(postdisp)
                    logdisp <- mixfit$meanprec
                    p0 <- mixfit$mixp[1]
                } else {
                    logdisp <- gaussMLwithmean(postdisp)
                }
            }
            diracprob <- c(p0,1-p0)
            
            rn <- seq(-8,8,by=0.01)
            rnprec <- seq(-8,1.5,by=0.01) #not
            if(!iselement2(inputpar$shrinkfixed,excludefornull)){
            KSfixed <- max(abs(pnorm(rn,mean=fixedprev[1],sd=sqrt(1/fixedprev[2]))-pnorm(rn,mean=fixed[1],sd=sqrt(1/fixed[2]))))
            } else {
            KSfixed <- max(c(
            abs(diracprobprev[2]*pnorm(rn[rn<0],mean=fixedprev[1],sd=sqrt(1/fixedprev[2]))-
            diracprob[2]*pnorm(rn[rn<0],mean=fixed[1],sd=sqrt(1/fixed[2]))),
            diracprobprev[2]*abs(pnorm(rn[rn>=0],mean=fixedprev[1],sd=sqrt(1/fixedprev[2])+diracprobprev[1])-(diracprob[2]*pnorm(rn[rn>=0],mean=fixed[1],sd=sqrt(1/fixed[2])) +
            +diracprob[1]))))    
            }
            whi <- which(!iselement2(inputpar$shrinkaddfixed,excludefornull))
            KSmaxaddfixed <- 0
            if(length(whi)>0){
                KSaddfixedall <- sapply(whi,function(i) {addfixedprevi <- addfixedprev[[i]];addfixedi <- addfixed[[i]];
                KSaddfixedi <- max(c(
                abs(diracprobprev[2]*pnorm(rn[rn<0],mean=addfixedprevi[1],sd=sqrt(1/addfixedprevi[2]))-
                diracprob[2]*pnorm(rn[rn<0],mean=addfixedi[1],sd=sqrt(1/addfixedi[2]))),
                diracprobprev[2]*abs(pnorm(rn[rn>=0],mean=addfixedprevi[1],sd=sqrt(1/addfixedprevi[2])+diracprobprev[1])-(diracprob[2]*pnorm(rn[rn>=0],mean=addfixedi[1],sd=sqrt(1/addfixedi[2])) +
                +diracprob[1])))) 
                })
            KSmaxaddfixed <- max(KSaddfixedall)
            }
            KSlogdisp <- max(abs(pnorm(rn,mean=logdispprev[1],sd=sqrt(1/logdispprev[2]))-pnorm(rn,mean=logdisp[1],sd=sqrt(1/logdisp[2]))))
            KSlogitp0 <- max(abs(pnorm(rn,mean=logitp0prev[1],sd=sqrt(1/logitp0prev[2]))-pnorm(rn,mean=logitp0[1],sd=sqrt(1/logitp0[2]))))
            
            if(!iselement2(inputpar$shrinkrandom,excludefornull) & !mixtrand){
            KSrandomprec <- max(abs(plgamma(rnprec,location=0,scale=1/randomprecprev[2],randomprecprev[1])-plgamma(rnprec,location=0,scale=1/randomprec[2],randomprec[1])))
            } else {
            KSrandomprec <- max(diracprobprev[2]*abs(plgamma(rnprec,location=0,scale=1/randomprecprev[2],randomprecprev[1])-
            diracprob[2]*plgamma(rnprec,location=0,scale=1/randomprec[2],randomprec[1]))) 
            }
            
            KS <- c(KSfixed,KSmaxaddfixed,KSlogdisp,KSlogitp0,KSrandomprec)
            names(KS) <- c("KSfixed","KSmaxaddfixed","KSlogdisp","KSlogitp0","KSrandomprec")
        ksall <- rbind(ksall,KS)
        KSmax <- max(KS[-length(KS)])
        #
        #create CURRENT list of relevant parameters
        pmlist <- list()
        fixed <- as.numeric(fixed)
        randomprec <- as.numeric(randomprec)
        logitp0 <- as.numeric(logitp0)
        logdisp <- as.numeric(logdisp)
        
        pmlist <- c(pmlist,list(mufixed=fixed[1],precfixed=fixed[2]))
        
        #NEW >= 1.8
        if(!is.null(shrinkaddfixed)){
        for(i in 1:els) {
            addf <- as.numeric(addfixed[[i]])
            name <- c(paste("mu",shrinkaddfixed[i],sep=""),
            paste("prec",shrinkaddfixed[i],sep="")) 
            toadd <- list(muaddfixed=addf[1],precaddfixed=addf[2])
            names(toadd) <- name
            pmlist <- c(pmlist,toadd)
        }
        } else {
        pmlist <- c(pmlist,list(muaddfixed=addfixed[[1]][1],precaddfixed=addfixed[[1]][2]))
        }
    
        pmlist <- c(pmlist,list(shaperand = randomprec[1],raterand=randomprec[2]))
        
        #NEW210
        if(!is.null(shrinkaddrandom)) {
        for(i in 1:elsr) {addr <- as.numeric(addrandom[[i]]);name <- c(paste("shape",shrinkaddrandom[i],sep=""),
        paste("rate",shrinkaddrandom[i],sep="")); toadd <- list(shapeaddr=addr[1],rateaddr=addr[2]);
        names(toadd) <- name;
        pmlist <- c(pmlist,toadd)}
        } else {
        pmlist <- c(pmlist,list(shapeaddr=addrandomprec[1],rateaddr=addrandomprec[2]))
        }
            
        pmlist <- c(pmlist,list(mup0=logitp0[1],precp0=logitp0[2])) 
    
        pmlist <- c(pmlist,list(mudisp=logdisp[1],precdisp=logdisp[2],mixp = diracprob))
    
        paranew <- c(unlist(pmlist),nfeat=ngenej, meanmlik = mlikmean)
        paraall <- rbind(paraall,paranew)
        }
        ################## MONITORING COVERGENCE ##########################
        print(paranew)
        print(paranew-paraprev)
        print(ksall)
        iter <- iter+1
        if(iter > maxiter | (KSmax <= tol & KSrandomprec<=tolrand) | mlikconv) moreiter <- FALSE  
        mlikprev <- mlikmean
        paraprev <- paranew 
        } #END OF WHILE LOOP
    } #END OF FOR LOOP 
#NEW >=1.8
#NEW210
    ret <- list(pmlist=pmlist,ksall=ksall,paraall=paraall,curvedispfun=curvedispfun,inputpar=inputpar,addfixed=addfixed,addrandom=addrandom,typelik="count")
return(ret)
} #END OF FUNCTION
