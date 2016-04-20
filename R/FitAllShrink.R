FitAllShrink <-
function(forms,dat,shrinksimul,finalprior=FALSE,dispersefixed=10, disperseaddfixed=1, disperserandom = 1, maxprecfixed=4, fams="zinb", ncpus = 
2, effoutput =TRUE, keepmargrand = FALSE, keepmarghyper=TRUE, setthreads1=TRUE,showupdate=FALSE,silentINLA=TRUE,updateby=5000,ndigits=5,
addpackage = NULL, safemode=TRUE,...){
#maxprecfixed: the maximum allowed precision in the prior for the main fixed parameter for fitting; useful when the 
#precision estimate of the simultaneous procedure is extremely high 
#forms = formTM;dat=mirnorm[whsig,];shrinksimul=shrinksimul

if(finalprior) {dispersefixed <- 1;disperseaddfixed<-1;disperserandom<-1}
if(dispersefixed != 1) print("IF YOU DO NOT INTEND TO USE EITHER OF THE FUNCTIONS MixtureUpdatePrior OR NonParaUpdatePrior 
PLEASE USE finalprior=TRUE")

typelik <- shrinksimul$typelik
ip <- shrinksimul$inputpar
curvedispfun <- shrinksimul$curvedispfun

if(typelik=="count"){
logdisp_shr <- c(mu=shrinksimul$pmlist$mudisp,prec=shrinksimul$pmlist$precdisp)
logitp0_shr <- c(mu=shrinksimul$pmlist$mup0,prec=shrinksimul$pmlist$precp0)
} else {
logdisp_shr <- c(0,0.01)
logitp0_shr <- c(0,0.01)
}

if(typelik=="gaussian"){
precerr_shr <- c(shapeerr=shrinksimul$pmlist$shapeerr,rateerr=shrinksimul$pmlist$rateerr)
} else {
precerr_shr <- c(0.001,0.001)
}

randeff_shr <- c(shape=shrinksimul$pmlist$shaperand/disperserandom,rate=shrinksimul$pmlist$raterand/disperserandom)

shrinkrandom <- ip$shrinkrandom
if(!is.null(shrinkrandom)) {
    if(mode(forms)=="call") form_shr <- randreplace(forms,shrinkrandom,randeff_shr) else 
    form_shr <- lapply(forms,randreplace,shrinkrandom=shrinkrandom,initrandomprec=randeff_shr)
    } else {
    form_shr <- forms
    }
    
shrinkaddrandom <- ip$shrinkaddrandom
addrandom <- shrinksimul$addrandom
    if(!is.null(shrinkaddrandom)){
    elsr <- length(shrinkaddrandom)
                for(i in 1:elsr){
                form_shr <- randreplace(form_shr,shrinkaddrandom[i],addrandom[[i]])
                }
                }
#IT IS WISE TO OVERDISPERSE THE PRIOR FOR FITTING IF ONE WANTS TO UPDATE THE PRIOR LATER

precfixeddisp <- shrinksimul$pmlist$precfixed/dispersefixed 
precfixeddisp <- min(maxprecfixed,precfixeddisp)
mufixed <- shrinksimul$pmlist$mufixed
shrinkfixed <- ip$shrinkfixed
shrinkaddfixed <- ip$shrinkaddfixed
addfixed <- shrinksimul$addfixed

addfixed <- lapply(addfixed,function(af) c(af[1],min(maxprecfixed,af[2]/disperseaddfixed)))
              
 cf <- list(mean=list(default=0), prec = list(default=0.01))
            if(!is.null(shrinkfixed)){
                if(is.factor(try(get(shrinkfixed),silent=TRUE))) shrinkfixed <- fact2vec(shrinkfixed)
                   nm <- names(cf[[1]])
                   for(j1 in 1:length(shrinkfixed)){       
                    cf$mean <- c(cf$mean,list(mufixed))
                    cf$prec <- c(cf$prec,list(precfixeddisp))
                    }
                    names(cf[[1]]) <- c(nm,shrinkfixed)
                    names(cf[[2]]) <- c(nm,shrinkfixed) 
                    } 

            if(!is.null(shrinkaddfixed)){ 
                for(i in 1:length(shrinkaddfixed)){
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

            
if((fams[1]=="poisson" | fams[1]=="zip") & (shrinksimul$pmlist$mixp[1] < 0.1)){
print("Estimated proportion to follow (ZI-)Poisson distribution is very low.")
print("You may want to abort this computation and continu with the results from the (ZI-)NB fit")
}
toret <- FitInlaAll(form_shr,dat, fams=fams, logdisp = logdisp_shr, precerr = precerr_shr, curvedispfun=curvedispfun, logitp0=logitp0_shr, 
ncpus = ncpus, effoutput = effoutput, keepmargrand = keepmargrand, keepmarghyper=keepmarghyper, setthreads1=setthreads1,
showupdate=showupdate, silentINLA=silentINLA, updateby=updateby,ndigits=ndigits,
addpackage = addpackage, cf=cf, safemode=safemode, ...)
priors <- list(mufixed=mufixed,precfixed=precfixeddisp,addfixed=addfixed,addrandom=addrandom,shaperand = randeff_shr[1],raterand=randeff_shr[2],
mudisp = logdisp_shr[1],precdisp = logdisp_shr[2],mup0 = logitp0_shr[1],precp0 = logitp0_shr[2], shapeerr = precerr_shr[1], rateerr = precerr_shr[2],
finalprior=finalprior)
return(list(res=toret,priors=priors))
}
