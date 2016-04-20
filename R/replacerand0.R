replacerand0 <-
function(form,shrinkrandom){
#shrinkrandom="tissue"; form <- y ~ 1 + cellfac + f(tissue) + f(tissue2) + f( seqno, model="z", Z=ZSpl(numKnots, deg, numtp, numcl, a, b), initial=3, param= randeff_init); 
frmchr <- as.character(form)[[3]]
sp <- strsplit(frmchr,"\\+")
sp <- as.vector(sapply(sp,function(tt){## remove whitespace
    gsub(" ","",tt)
    }))
sr1 <- paste("\\(",shrinkrandom,",",sep="")  
sr2 <- paste("\\(",shrinkrandom,"\\)",sep="")   
wh <- union(union(which(sp==shrinkrandom),grep(sr1,sp)),grep(sr2,sp))
if(is.numeric(wh)) {
    sp <- sp[-wh]
    fnew <- sp[1]; if(length(sp)>1) {for(i in 2:length(sp)) fnew <- paste(fnew,sp[i],sep="+")}
    formret <- formula(paste("y ~ ",fnew,sep=""))
return(formret)
} else {
return(form)
}
}
