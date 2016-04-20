ScreenData <-
function(dat,variable,np=TRUE,ncpus=2,updateby=10000){
el <- nrow(dat)
pv <- rep(NA,el)
if(el>updateby) sec <- seq(1,el,by=updateby) else sec <- c(1,el)
secl <- length(sec)
if(sec[secl] > el-2) sec[secl] <- el+1 else sec <- c(sec,el+1)
secl <- length(sec)
if(ncpus>1) sfInit(parallel=TRUE,cpus=ncpus) 
for(k in 1:(secl-1)){
thedatak <- dat[(sec[k]:(sec[k+1]-1)),,drop=FALSE]
pvk <- screendatahelp(thedatak,variable=variable,np=np,ncpus=ncpus)
pv[(sec[k]):(sec[k+1]-1)] <- pvk
print(paste((sec[k+1]-1),"Cases done (",round(100*(sec[k+1]-1)/el),"%)"))
}
if(ncpus>1) sfStop()
return(pv)
}
